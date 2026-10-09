#----------------------------------------------------------#
# 04 Publication figures -----
#----------------------------------------------------------#
# Input:  Data/Processed/fields_clean.rds, Outputs/Models/spatial.rds
# Output: Outputs/Figures/*.{pdf,png,tiff}
#
# Effects are drawn from the final model of each hypothesis (spatial refit
# where residual autocorrelation required one), at population level (random
# effects at zero), with 95% CIs built on the link scale so they cannot cross
# zero. Count models are shown as density per 1 km of transect; presence
# models as detection probability on a 5-km transect.

source(here::here("R", "00_setup.R"))

fields <- readRDS(path_clean_rds)
spatial <- readRDS(path_spatial_rds)

final <- spatial %>%
  filter(role == "main") %>%
  select(model_id, response, fit_final, model_data, sim_final, use_spatial)

get_model <- function(id) final %>% filter(model_id == id)

density_label <- paste0(str_to_sentence(count_unit), " per km of transect")
presence_label <- paste0("Detection probability (", line_reference_km, "-km transect)")

#----------------------------------------------------------#
# Prediction helper -----
#----------------------------------------------------------#
# Other continuous predictors are held at their sample mean and time of day
# at "morning" (the most common category).
predict_effect <- function(id, term, n = 200) {
  m <- get_model(id)
  fit <- m$fit_final[[1]]
  dat <- m$model_data[[1]]
  info <- predictor_info %>% filter(term == !!term)
  stopifnot(length(fixef(fit)$zi) == 0)

  raw_seq <- seq(min(dat[[info$raw]]), max(dat[[info$raw]]), length.out = n)

  new_data <- add_spatial_coords(dat)[rep(1, n), ]
  for (t in intersect(predictor_info$term, names(dat))) {
    new_data[[t]] <- mean(dat[[t]], na.rm = TRUE)
  }
  if ("time_of_day_2" %in% names(new_data)) {
    new_data$time_of_day_2 <- factor("morning", levels = levels(dat$time_of_day_2))
  }
  new_data[[term]] <- (raw_seq - info$centre) / info$scale
  new_data$x <- raw_seq

  is_presence <- m$response == "presence"
  new_data$line_km <- if (is_presence) line_reference_km else 1

  p <- predict(fit, newdata = new_data, re.form = NA, se.fit = TRUE, type = "link")
  inv <- if (is_presence) \(eta) 1 - exp(-exp(eta)) else exp

  tibble(
    x = raw_seq,
    fit = inv(p$fit),
    low = inv(p$fit - qnorm(0.975) * p$se.fit),
    high = inv(p$fit + qnorm(0.975) * p$se.fit)
  )
}

observed_points <- function(id, term) {
  m <- get_model(id)
  dat <- m$model_data[[1]]
  raw <- predictor_info$raw[predictor_info$term == term]
  y <- if (m$response == "presence") dat$presence else dat[[m$response]] / dat$line_km
  tibble(x = dat[[raw]], y = y)
}

#----------------------------------------------------------#
# Panel builder -----
#----------------------------------------------------------#
sqrt_breaks <- c(0, 0.1, 0.25, 0.5, 1, 1.5, 2)
# Built from month.abb rather than format(), which follows the system locale
date_labels <- \(x) {
  d <- as_date(x - 1, origin = "2025-01-01")
  if_else(is.na(d), NA_character_, paste(mday(d), month.abb[month(d)]))
}
date_breaks <- \(lims) {
  days <- seq(as_date("2025-01-01"), as_date("2025-12-31"), by = "day")
  # Mid-month ticks only where the axis is short enough to keep labels apart
  mid <- if (diff(lims) > 70) 1 else c(1, 15)
  doy <- yday(days[mday(days) %in% mid])
  doy[doy >= lims[1] & doy <= lims[2]]
}

effect_panel <- function(id, term, x_jitter = 0) {
  m <- get_model(id)
  pred <- predict_effect(id, term)
  obs <- observed_points(id, term)
  info <- predictor_info %>% filter(term == !!term)
  is_presence <- m$response == "presence"

  p <- ggplot(pred, aes(x = x)) +
    geom_point(
      data = obs, aes(x = x, y = y),
      position = position_jitter(width = x_jitter, height = if (is_presence) 0.03 else 0, seed = seed),
      colour = col_points, alpha = 0.3, size = 0.7, stroke = 0
    ) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = col_main, alpha = 0.25) +
    geom_line(aes(y = fit), colour = col_main, linewidth = 0.7) +
    labs(x = info$axis_label, y = if (is_presence) presence_label else density_label)

  p <- if (is_presence) {
    p + scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.25))
  } else {
    # Square-root axis: most fields have zero or low density, so a linear axis
    # would compress the fitted line against the floor.
    p + scale_y_continuous(trans = "sqrt", breaks = sqrt_breaks)
  }
  if (str_detect(term, "^doy")) p <- p + scale_x_continuous(breaks = date_breaks, labels = date_labels)
  p
}

time_of_day_panel <- function(id) {
  m <- get_model(id)
  fit <- m$fit_final[[1]]
  dat <- m$model_data[[1]]
  levs <- levels(droplevels(dat$time_of_day_2))
  new_data <- add_spatial_coords(dat)[rep(1, length(levs)), ]
  for (t in intersect(predictor_info$term, names(dat))) {
    new_data[[t]] <- mean(dat[[t]], na.rm = TRUE)
  }
  new_data$time_of_day_2 <- factor(levs, levels = levels(dat$time_of_day_2))
  new_data$line_km <- 1
  p <- predict(fit, newdata = new_data, re.form = NA, se.fit = TRUE, type = "link")
  est <- tibble(
    level = factor(levs, levels = levs),
    n = as.integer(table(droplevels(dat$time_of_day_2))[levs]),
    fit = exp(p$fit),
    low = exp(p$fit - qnorm(0.975) * p$se.fit),
    high = exp(p$fit + qnorm(0.975) * p$se.fit)
  )
  ggplot(est, aes(x = level, y = fit)) +
    geom_pointrange(aes(ymin = low, ymax = high), colour = col_main, size = 0.3, linewidth = 0.6) +
    geom_text(aes(y = high, label = paste0("n = ", n)), vjust = -0.6, size = 2.3, colour = "grey30") +
    scale_x_discrete(labels = \(x) str_replace_all(x, c("/" = "/\n", " \\+ " = " +\n"))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = "Time of day, visit 2", y = density_label)
}

tag_panels <- function(p) p + plot_annotation(tag_levels = "A")

#----------------------------------------------------------#
# Figures 1-5: hypothesis effects -----
#----------------------------------------------------------#
fig_1 <- tag_panels(
  effect_panel("H1", "spruce_s") +
    effect_panel("H2", "ndop_s", x_jitter = 0.15) +
    effect_panel("H3", "prob_s")
)
save_figure(fig_1, "fig1_h1_h3_habitat", height_mm = 65)

fig_2 <- tag_panels(
  effect_panel("H4", "calamity_s") +
    effect_panel("H5", "dry_s") +
    effect_panel("H6", "harvest_s")
)
save_figure(fig_2, "fig2_h4_h6_disturbance", height_mm = 65)

fig_3 <- tag_panels(
  effect_panel("H7", "doy_1_s", x_jitter = 0.3) +
    effect_panel("H8", "doy_2_s", x_jitter = 0.3) +
    time_of_day_panel("H8")
)
save_figure(fig_3, "fig3_h7_h8_survey_date", height_mm = 70)

fig_4 <- tag_panels(
  (effect_panel("H9_v1", "loudness_s", x_jitter = 0.3) + labs(subtitle = "Visit 1, provoked")) +
    (effect_panel("H9_v2", "loudness_s", x_jitter = 0.3) + labs(subtitle = "Visit 2, provoked"))
)
save_figure(fig_4, "fig4_h9_device_loudness", width_mm = width_double, height_mm = 70)

fig_5 <- tag_panels(
  (effect_panel("H10_H11_pres", "lat_s") + effect_panel("H10_H11_pres", "lon_s")) /
    (effect_panel("H10_H11_dens", "lat_s") + effect_panel("H10_H11_dens", "lon_s"))
)
save_figure(fig_5, "fig5_h10_h11_geography", height_mm = 120)

#----------------------------------------------------------#
# Maps -----
#----------------------------------------------------------#
# Bundled 1:50m Natural Earth borders, so the map needs no download
czechia <- st_as_sf(rnaturalearthdata::countries50) %>%
  filter(adm0_a3 == "CZE") %>%
  st_transform(crs_sjtsk)

# Each field is a KFME sixteenth-cell: 2.5' of longitude by 1.5' of latitude
cell_polygon <- function(lon, lat) {
  dx <- 1 / 48
  dy <- 0.0125
  st_polygon(list(rbind(
    c(lon - dx, lat - dy), c(lon + dx, lat - dy), c(lon + dx, lat + dy),
    c(lon - dx, lat + dy), c(lon - dx, lat - dy)
  )))
}

# Added last: geom_sf() sets its own default coordinate system
map_coord <- coord_sf(crs = crs_sjtsk, datum = crs_wgs84, expand = FALSE)

field_cells <- fields %>%
  mutate(geometry = st_sfc(map2(lon, lat, cell_polygon), crs = crs_wgs84)) %>%
  st_as_sf() %>%
  st_transform(crs_sjtsk)

map_base <- function() {
  list(
    geom_sf(data = czechia, fill = "grey97", colour = "grey40", linewidth = 0.3),
    theme(
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
      axis.line = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",
      legend.title.position = "top",
      legend.title = element_text(hjust = 0.5),
      legend.key.width = unit(12, "mm"),
      legend.key.height = unit(2.5, "mm")
    )
  )
}

fig_6 <- ggplot() +
  map_base() +
  geom_sf(data = field_cells, aes(fill = density), colour = "grey30", linewidth = 0.05) +
  scale_fill_viridis_c(
    name = density_label, option = "viridis", direction = -1,
    trans = "sqrt", breaks = sqrt_breaks, na.value = "grey70"
  ) +
  map_coord
save_figure(fig_6, "fig6_map_density", height_mm = 105)

residual_cells <- function(id) {
  m <- get_model(id)
  m$model_data[[1]] %>%
    transmute(kfme_cell, residual = residuals(m$sim_final[[1]])) %>%
    inner_join(field_cells %>% select(kfme_cell), by = "kfme_cell") %>%
    st_as_sf()
}

residual_map <- function(id, title) {
  ggplot() +
    map_base() +
    geom_sf(data = residual_cells(id), aes(fill = residual), colour = "grey30", linewidth = 0.05) +
    scale_fill_distiller(
      name = "Scaled residual (DHARMa)", palette = "PuOr", direction = 1,
      limits = c(0, 1), breaks = c(0, 0.5, 1)
    ) +
    map_coord +
    labs(subtitle = title)
}

fig_7 <- tag_panels(
  (residual_map("H10_H11_pres", "Presence model (latitude + longitude)") +
     residual_map("H10_H11_dens", "Density model (latitude + longitude)")) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
)
save_figure(fig_7, "fig7_map_residuals", height_mm = 80)

#----------------------------------------------------------#
# Supplementary: residual correlograms for every main model -----
#----------------------------------------------------------#
correlograms <- spatial %>%
  filter(role == "main") %>%
  select(model_id, correlogram_before, correlogram_after) %>%
  pivot_longer(-model_id, names_to = "fit_type", values_to = "corr") %>%
  filter(!map_lgl(corr, is.null)) %>%
  unnest(corr) %>%
  mutate(
    fit_type = recode(fit_type, correlogram_before = "Non-spatial", correlogram_after = "Spatial (exponential)"),
    significant = p_value < 0.05,
    model_id = factor(model_id, levels = unique(spatial$model_id))
  )

fig_s1 <- ggplot(correlograms, aes(distance_km, moran_I, colour = fit_type)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_point(aes(shape = significant), size = 1.4, fill = "white") +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 19), labels = c("p ≥ 0.05", "p < 0.05"), name = NULL) +
  scale_colour_manual(values = c(col_main, col_accent), name = NULL) +
  facet_wrap(~model_id, ncol = 4) +
  scale_x_continuous(breaks = c(0, 50, 100)) +
  labs(x = "Distance class midpoint (km)", y = "Moran's I of residuals") +
  theme(legend.position = "bottom")
save_figure(fig_s1, "figS1_residual_correlograms", height_mm = 190)

message("04: figures written to ", dir_figures)
