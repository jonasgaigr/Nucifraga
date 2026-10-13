#----------------------------------------------------------#
# 06 Exploratory PCA of the habitat variables -----
#----------------------------------------------------------#
# Input:  Data/Processed/fields_clean.rds
# Output: Outputs/Tables/pca_{variance,loadings,density_model}.csv,
#         Outputs/Figures/figS3_pca_habitat.{pdf,png,tiff}
#
# Exploratory only; the hypothesis tests are in 02-03. Replaces the earlier
# vegan PCA, which mixed responses (total and spontaneous counts) with
# predictors, used raw m² areas, included calamity together with its own
# components, and grouped points by field ID.
#
# Here the PCA describes how the habitat predictors co-vary. Nutcracker density
# and the coordinates are passive: they are projected onto the ordination but
# do not shape it. A density model on the first two components, with the same
# effort offset (and spatial refit) as the main models, shows how much of the
# habitat signal the joint gradients carry.

source(here::here("R", "00_setup.R"))

fields <- readRDS(path_clean_rds)

# Same fields as the total-count models, so density and the PC model are
# comparable with H1-H6
pca_data <- fields %>% filter(!is.na(line_km), visit_2_done)

#----------------------------------------------------------#
# Active and passive variables -----
#----------------------------------------------------------#
# Calamity area is exactly dry + harvested area, so it is left out: it would
# only duplicate its components. Transformations remove the strong right skew
# (NDOP skewness 4.3, occurrence probability 2.2, dry 5.1, harvest 1.9);
# all proportions are strictly between 0 and 1.
logit <- \(p) log(p / (1 - p))

active <- tribble(
  ~variable,      ~label,                     ~transform,
  "spruce_km2",   "Spruce forest",            identity,
  "ndop_records", "NDOP records",             log1p,
  "prob_occ",     "Occurrence probability",   log,
  "prop_dry",     "Standing dead spruce",     logit,
  "prop_harvest", "Salvage-logged",           logit
)

active_matrix <- map2(active$variable, active$transform, \(v, f) f(pca_data[[v]])) %>%
  set_names(active$label) %>%
  as_tibble() %>%
  as.matrix()
stopifnot(all(is.finite(active_matrix)))

passive <- tibble(
  `Density (√)` = sqrt(pca_data$density),
  Latitude = pca_data$lat,
  Longitude = pca_data$lon
)

#----------------------------------------------------------#
# PCA -----
#----------------------------------------------------------#
# Correlation-matrix PCA: the variables are on different scales
pca <- prcomp(active_matrix, center = TRUE, scale. = TRUE)
n_pc <- ncol(active_matrix)

# Component signs are arbitrary. PC1 is oriented so that more spruce forest is
# positive; every other component so that its largest loading is positive.
orientation <- sign(pca$rotation[cbind(apply(abs(pca$rotation), 2, which.max), seq_len(n_pc))])
orientation[1] <- sign(pca$rotation["Spruce forest", 1])
pca$rotation <- sweep(pca$rotation, 2, orientation, `*`)
pca$x <- sweep(pca$x, 2, orientation, `*`)

# Broken-stick expectation: components whose share exceeds it carry more
# structure than random partitioning of the variance would. Components are
# retained up to the first one that falls below it.
broken_stick <- map_dbl(seq_len(n_pc), \(k) sum(1 / (k:n_pc)) / n_pc)

pca_variance <- tibble(
  component = paste0("PC", seq_len(n_pc)),
  eigenvalue = pca$sdev^2,
  prop_variance = eigenvalue / sum(eigenvalue),
  cum_variance = cumsum(prop_variance),
  broken_stick = broken_stick,
  retained = cumall(prop_variance > broken_stick)
)

# Correlation-biplot coordinates: correlation of each variable with the
# standardised scores (for active variables equal to loading × sdev)
scores <- pca$x %*% diag(1 / pca$sdev)
colnames(scores) <- colnames(pca$x)

loadings <- bind_rows(
  as_tibble(cor(active_matrix, scores), rownames = "variable") %>% mutate(role = "active"),
  as_tibble(cor(passive, scores), rownames = "variable") %>% mutate(role = "passive")
) %>%
  select(role, variable, PC1:PC3)

write_csv(pca_variance, file.path(dir_tables, "pca_variance.csv"))
write_csv(loadings, file.path(dir_tables, "pca_loadings.csv"))

#----------------------------------------------------------#
# Density on the habitat gradients -----
#----------------------------------------------------------#
model_data <- pca_data %>%
  mutate(pc1 = scores[, "PC1"], pc2 = scores[, "PC2"])

fit_nb <- glmmTMB(build_formula("count_total", c("pc1", "pc2")), family = nbinom2, data = model_data)
fit_ns <- if (converged(fit_nb) && sigma(fit_nb) < 1e3) fit_nb else update(fit_nb, family = poisson)

# The main models all showed residual autocorrelation, so the spatial refit is
# reported alongside
fit_sp <- fit_spatial(fit_ns, model_data, "count_total", c("pc1", "pc2"))

pc_model <- bind_rows(
  tidy(fit_ns, effects = "fixed", conf.int = TRUE) %>% mutate(fit_type = "non-spatial"),
  tidy(fit_sp, effects = "fixed", conf.int = TRUE) %>% mutate(fit_type = "spatial")
) %>%
  filter(component == "cond", term != "(Intercept)") %>%
  transmute(
    fit_type, family = c(family(fit_ns)$family, family(fit_sp)$family)[match(fit_type, c("non-spatial", "spatial"))],
    term, step = "per 1 SD of component score", n = nrow(model_data),
    estimate_log = estimate, se = std.error,
    irr = exp(estimate), irr_ci_low = exp(conf.low), irr_ci_high = exp(conf.high),
    p_value = p.value, converged = c(converged(fit_ns), converged(fit_sp))[match(fit_type, c("non-spatial", "spatial"))]
  )
write_csv(pc_model, file.path(dir_tables, "pca_density_model.csv"))

#----------------------------------------------------------#
# Figure S3 -----
#----------------------------------------------------------#
pc_label <- \(k) paste0("PC", k, " (", round(100 * pca_variance$prop_variance[k]), " %)")

# Correlation circle: arrows are correlations with the scores, so their true
# lengths are comparable. Latitude and longitude are weakly related to PC1-2
# (|r| < 0.2) and stay in pca_loadings.csv rather than on the plot.
short_label <- c(`Occurrence probability` = "Occurrence prob.", `Standing dead spruce` = "Standing dead")

# Pushes labels on the same side of the plot apart vertically, so arrows that
# point in similar directions keep readable labels
spread_labels <- function(y, min_gap) {
  o <- order(y)
  ys <- y[o]
  for (i in seq_along(ys)[-1]) ys[i] <- max(ys[i], ys[i - 1] + min_gap)
  y[o] <- ys
  y
}

circle_arrows <- loadings %>%
  filter(variable %in% c(active$label, "Density (√)")) %>%
  mutate(
    label = coalesce(short_label[variable], variable),
    lx = PC1 * 1.08 + if_else(PC1 >= 0, 0.03, -0.03),
    ly = PC2 * 1.08,
    hjust = if_else(PC1 >= 0, 0, 1)
  ) %>%
  group_by(hjust) %>%
  mutate(ly = spread_labels(ly, 0.12)) %>%
  ungroup()

unit_circle <- tibble(a = seq(0, 2 * pi, length.out = 200), x = cos(a), y = sin(a))

axis_lines <- list(
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3),
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3)
)

ordination_theme <- theme(
  panel.grid.major.y = element_blank(),
  legend.position = "bottom",
  legend.title.position = "top",
  legend.title = element_text(hjust = 0.5)
)

panel_scree <- ggplot(pca_variance, aes(x = component)) +
  geom_col(aes(y = prop_variance), fill = col_main, width = 0.7) +
  geom_point(aes(y = broken_stick), shape = 21, fill = "white", colour = "grey15", size = 1.5) +
  scale_y_continuous(labels = \(x) paste0(round(100 * x), " %"), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Variance explained")

panel_circle <- ggplot(circle_arrows) +
  axis_lines +
  geom_path(data = unit_circle, aes(x, y), colour = "grey70", linewidth = 0.3) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2, colour = role, linetype = role),
               arrow = arrow(length = unit(1.2, "mm")), linewidth = 0.45) +
  geom_text(aes(x = lx, y = ly, label = label, hjust = hjust, colour = role), size = 2.3) +
  scale_colour_manual(values = c(active = "grey15", passive = col_accent), guide = "none") +
  scale_linetype_manual(values = c(active = "solid", passive = "22"), guide = "none") +
  coord_equal(xlim = c(-1.75, 1.75), ylim = c(-1.1, 1.1)) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  labs(x = pc_label(1), y = pc_label(2)) +
  ordination_theme

score_points <- model_data %>%
  transmute(PC1 = pc1, PC2 = pc2, density, selected_2025)
# Shared limits for C and D that also contain the selection ellipses
ellipses <- bind_rows(map(c("no", "yes"), \(lev) {
  d <- score_points %>% filter(selected_2025 == lev) %>% select(PC1, PC2)
  ctr <- colMeans(d)
  shape <- chol(cov(d))
  radius <- sqrt(2 * qf(0.95, 2, nrow(d) - 1))
  a <- seq(0, 2 * pi, length.out = 100)
  as_tibble(t(ctr + t(radius * cbind(cos(a), sin(a)) %*% shape)))
}))
score_xlim <- range(score_points$PC1, ellipses$PC1)
score_ylim <- range(score_points$PC2, ellipses$PC2)

panel_density <- ggplot(score_points %>% arrange(density), aes(PC1, PC2)) +
  axis_lines +
  geom_point(aes(fill = density), shape = 21, colour = "grey40", stroke = 0.15, size = 1.2) +
  scale_fill_gradient(name = paste0(str_to_sentence(count_unit), " per km of transect"),
                      low = pal_sequential[1], high = pal_sequential[2], trans = "sqrt",
                      breaks = c(0, 0.25, 0.5, 1, 1.5)) +
  coord_equal(xlim = score_xlim, ylim = score_ylim) +
  labs(x = pc_label(1), y = pc_label(2)) +
  ordination_theme +
  theme(legend.key.width = unit(10, "mm"), legend.key.height = unit(2.5, "mm"))

panel_selection <- ggplot(score_points, aes(PC1, PC2, colour = selected_2025, shape = selected_2025)) +
  axis_lines +
  geom_point(size = 1, alpha = 0.6, stroke = 0) +
  stat_ellipse(aes(group = selected_2025), level = 0.95, linewidth = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = c(no = col_main, yes = col_accent), name = "Selected for 2025",
                      labels = c(no = "No", yes = "Yes")) +
  scale_shape_manual(values = c(no = 16, yes = 17), name = "Selected for 2025",
                     labels = c(no = "No", yes = "Yes")) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  coord_equal(xlim = score_xlim, ylim = score_ylim) +
  labs(x = pc_label(1), y = pc_label(2)) +
  ordination_theme

fig_s3 <- (panel_scree + panel_circle) / (panel_density + panel_selection) +
  plot_layout(heights = c(0.8, 1)) +
  plot_annotation(tag_levels = "A")
save_figure(fig_s3, "figS3_pca_habitat", width_mm = width_double, height_mm = 165)

message(
  "06: PCA on ", nrow(model_data), " fields; PC1 ", round(100 * pca_variance$prop_variance[1]),
  " %, PC2 ", round(100 * pca_variance$prop_variance[2]), " % of variance; ",
  sum(pca_variance$retained), " component(s) above the broken stick."
)
