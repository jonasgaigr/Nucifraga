data <- read_csv2(
  "Data/Processed/data_nuc.csv"
)
glimpse(data)
# Response matrix: numeric soil properties
vars <- data %>%
  dplyr::select(
    Celkem_zjisteno_oresniku,
    Pocet_nalezu,
    SMRK_m2,
    prob_1h_mean,
    Plocha_kurovcove_kalamity_celkem_m2,
    Plocha_souse_m2,
    Plocha_tezby_m2,
    Spontanne_1
)

# Verify no infinite or NA values remain
stopifnot(!any(is.na(vars)), !any(is.infinite(as.matrix(vars))))

# --- PCA (unconstrained ordination) ------------------------------------------
pca <- rda(vars, scale = TRUE)

# Fit abundance variables as vectors
ef_pca <- envfit(pca, vars, permutations = 999)

# Summary of PCA
summary(pca)

# Screeplot of eigenvalues
screeplot(pca, bstick = TRUE, main = "PCA – Nucifraga abundance")

# --- PCA site scores ---------------------------------------------------------
site_scores_df <- as.data.frame(scores(pca, display = "sites", scaling = 2))
site_scores_df$sample_place <- data$ID

# Arrows (biplot scores for abundance variables)
arrow_scores <- as.data.frame(scores(ef_pca, display = "vectors", scaling = 2))

# Compute centroids by site (if applicable)
centroids <- aggregate(cbind(PC1 = site_scores_df$PC1, PC2 = site_scores_df$PC2),
                       by = list(sample_place = site_scores_df$sample_place),
                       FUN = mean)

# --- PCA plot ---------------------------------------------------------------
pca_plot <- ggplot(site_scores_df, aes(x = PC1, y = PC2, color = sample_place)) +
  geom_point(size = 2) +
  geom_mark_hull(
    aes(group = sample_place, fill = sample_place),
    alpha = 0.2, concavity = 5, expand = unit(2, "mm"), show.legend = FALSE
  ) +
  geom_segment(data = arrow_scores,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               inherit.aes = FALSE,
               color = "red") +
  geom_text_repel(data = centroids,
                  aes(x = PC1, y = PC2, label = sample_place),
                  size = 4, fontface = "bold", color = "black") +
  geom_text_repel(data = arrow_scores,
                  aes(x = PC1, y = PC2, label = rownames(arrow_scores)),
                  inherit.aes = FALSE,
                  color = "red", size = 3) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("PCA (correlation biplot) – Nucifraga abundance")

ggsave("Outputs/Plots/pca_nucifraga.png", plot = pca_plot, height = 5, width = 6)
ggsave("Outputs/Plots/pca_nucifraga.pdf", plot = pca_plot, height = 5, width = 6, device = cairo_pdf)

# --- RDA (constrained ordination) -------------------------------------------
# Use explanatory variable (e.g., habitat type or site)
if ("ID" %in% names(data)) {
  rda_model <- rda(vars ~ sample_place, data = data, scale = TRUE)
  
  summary(rda_model)
  anova(rda_model)
  anova(rda_model, by = "axis")
  anova(rda_model, by = "terms")
  
  # Site scores
  site_scores_df <- as.data.frame(scores(rda_model, display = "sites", scaling = 2))
  site_scores_df$sample_place <- data$ID
  
  # Arrows (biplot scores)
  arrow_scores_rda <- as.data.frame(scores(rda_model, display = "bp", scaling = 2))
  arrow_scores_rda$varname <- rownames(arrow_scores_rda)
  colnames(arrow_scores_rda)[1:2] <- c("RDA1", "RDA2")
  
  # Centroids
  centroids <- aggregate(cbind(RDA1 = site_scores_df$RDA1, RDA2 = site_scores_df$RDA2),
                         by = list(sample_place = site_scores_df$sample_place),
                         FUN = mean)
  
  # Plot
  rda_plot <- ggplot(site_scores_df, aes(x = RDA1, y = RDA2, color = sample_place)) +
    geom_point(size = 2) +
    stat_chull(aes(group = sample_place), alpha = 0.2, geom = "polygon") +
    geom_text_repel(data = centroids,
                    aes(x = RDA1, y = RDA2, label = ID),
                    size = 4, fontface = "bold", color = "black") +
    theme_minimal() +
    ggtitle("RDA – Nucifraga abundance ~ sample_place") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "RDA1", y = "RDA2")
  
  ggsave("Outputs/Plots/rda_nucifraga.png", plot = rda_plot, height = 5, width = 6)
  ggsave("Outputs/Plots/rda_nucifraga.pdf", plot = rda_plot, height = 5, width = 6, device = cairo_pdf)
}
