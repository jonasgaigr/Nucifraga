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
#
# Figures 1-5 (and their subset versions for the poster) are written in
# English and in Czech (file names ending in _cs).

source(here::here("R", "00_setup.R"))

fields <- readRDS(path_clean_rds)
spatial <- readRDS(path_spatial_rds)

final <- spatial %>%
  select(model_id, response, fit_final, model_data, sim_final, use_spatial)

get_model <- function(id) final %>% filter(model_id == id)

#----------------------------------------------------------#
# Labels in English and Czech -----
#----------------------------------------------------------#
axis_labels <- tribble(
  ~term,        ~en,                                                                  ~cs,
  "spruce_s",   "Spruce stand area (km²)",                                            "Plocha smrkového porostu (km²)",
  "ndop_s",     "Number of records in NDOP (2014–2024)",                              "Počet záznamů v NDOP (2014–2024)",
  "prob_s",     "Probability of detection within one hour (2014–2017)",               "Pravděpodobnost zastižení během jedné hodiny (2014–2017)",
  "calamity_s", "Proportion of disturbed area to total spruce stand area",            "Podíl kalamitní plochy na celkovou plochu smrkového porostu",
  "dry_s",      "Proportion of unsalvaged disturbed area to total spruce stand area", "Podíl nevytěžené kalamitní plochy na celkovou plochu smrkového porostu",
  "harvest_s",  "Proportion of salvage-logged disturbed area to total spruce stand area", "Podíl vytěžené kalamitní plochy na celkovou plochu smrkového porostu",
  "doy_1_s",    "Date of 1st survey visit",                                           "Datum 1. kontroly",
  "doy_2_s",    "Date of 2nd survey visit",                                           "Datum 2. kontroly",
  "loudness_s", "Loudness of the playback device (dB at 1 m)",                        "Hlasitost zařízení určeného k provokaci (dB v 1 m)",
  "lat_s",      "Latitude (°N)",                                                      "Zeměpisná šířka (° s. š.)",
  "lon_s",      "Longitude (°E)",                                                     "Zeměpisná délka (° v. d.)"
)

text_labels <- list(
  en = list(
    density = "Number of territories per 1 km",
    presence = paste0("Detection probability (", line_reference_km, "-km transect)"),
    time_of_day = "Time of day, 2nd survey visit",
    time_of_day_levels = c(morning = "Morning", `afternoon/evening` = "Afternoon/evening",
                           `morning + evening` = "Morning + evening"),
    loudness_levels = c(`<80` = "< 80", `80-90` = "80–90", `90-100` = "90–100", `>=100` = "≥ 100"),
    visit_1 = "1st survey visit",
    visit_2 = "2nd survey visit",
    spontaneous = "Nutcrackers detected without playback",
    provoked = "Nutcrackers detected with playback",
    selected = "Fields selected for 2025 (Vyber_2025)",
    selected_low_calamity = "Fields selected for 2025, calamity area < 50% of spruce area",
    date = \(day, month) paste(day, month.abb[month])
  ),
  cs = list(
    density = "Počet teritorií / 1 km",
    presence = paste0("Pravděpodobnost zjištění (transekt ", line_reference_km, " km)"),
    time_of_day = "Denní doba, 2. kontrola",
    time_of_day_levels = c(morning = "Ráno", `afternoon/evening` = "Odpoledne/večer",
                           `morning + evening` = "Ráno + večer"),
    loudness_levels = c(`<80` = "< 80", `80-90` = "80–90", `90-100` = "90–100", `>=100` = "≥ 100"),
    visit_1 = "1. kontrola",
    visit_2 = "2. kontrola",
    spontaneous = "Ořešníci zjištění bez provokace",
    provoked = "Ořešníci zjištění provokací",
    selected = "Vylosovaná pole (Vyber_2025)",
    selected_low_calamity = "Vylosovaná pole, kalamita < 50 % plochy smrku",
    date = \(day, month) paste0(day, ". ", month, ".")
  )
)

# English labels for the maps and supplementary figures
density_label <- paste0(str_to_sentence(count_unit), " per km of transect")

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
# Panel builders -----
#----------------------------------------------------------#
sqrt_breaks <- c(0, 0.1, 0.25, 0.5, 1, 1.5, 2)
# Built from the label table rather than format(), which follows the system
# locale
date_labels <- \(lang) \(x) {
  d <- as_date(x - 1, origin = "2025-01-01")
  if_else(is.na(d), NA_character_, text_labels[[lang]]$date(mday(d), month(d)))
}
date_breaks <- \(lims) {
  days <- seq(as_date("2025-01-01"), as_date("2025-12-31"), by = "day")
  # Mid-month ticks only where the axis is short enough to keep labels apart
  mid <- if (diff(lims) > 70) 1 else c(1, 15)
  doy <- yday(days[mday(days) %in% mid])
  doy[doy >= lims[1] & doy <= lims[2]]
}

# The disturbance labels are long; wrapped to fit one third of a page width
axis_label <- function(term, lang) str_wrap(axis_labels[[lang]][axis_labels$term == term], 42)

effect_panel <- function(id, term, lang = "en", x_jitter = 0) {
  m <- get_model(id)
  pred <- predict_effect(id, term)
  obs <- observed_points(id, term)
  lab <- text_labels[[lang]]
  is_presence <- m$response == "presence"

  p <- ggplot(pred, aes(x = x)) +
    geom_point(
      data = obs, aes(x = x, y = y),
      position = position_jitter(width = x_jitter, height = if (is_presence) 0.03 else 0, seed = seed),
      colour = col_points, alpha = 0.3, size = 0.7, stroke = 0
    ) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = col_main, alpha = 0.25) +
    geom_line(aes(y = fit), colour = col_main, linewidth = 0.7) +
    labs(x = axis_label(term, lang), y = if (is_presence) lab$presence else lab$density)

  p <- if (is_presence) {
    p + scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.25))
  } else {
    # Square-root axis: most fields have zero or low density, so a linear axis
    # would compress the fitted line against the floor.
    p + scale_y_continuous(trans = "sqrt", breaks = sqrt_breaks)
  }
  if (str_detect(term, "^doy")) p <- p + scale_x_continuous(breaks = date_breaks, labels = date_labels(lang))
  p
}

# Predicted density per level of a factor, other predictors at their mean
# and time of day at morning
factor_panel <- function(id, variable, level_labels, x_label, lang = "en") {
  m <- get_model(id)
  fit <- m$fit_final[[1]]
  dat <- m$model_data[[1]]
  levs <- levels(droplevels(dat[[variable]]))
  new_data <- add_spatial_coords(dat)[rep(1, length(levs)), ]
  for (t in intersect(predictor_info$term, names(dat))) {
    new_data[[t]] <- mean(dat[[t]], na.rm = TRUE)
  }
  if (variable != "time_of_day_2" && "time_of_day_2" %in% names(new_data)) {
    new_data$time_of_day_2 <- factor("morning", levels = levels(dat$time_of_day_2))
  }
  new_data[[variable]] <- factor(levs, levels = levels(dat[[variable]]))
  new_data$line_km <- 1
  p <- predict(fit, newdata = new_data, re.form = NA, se.fit = TRUE, type = "link")
  est <- tibble(
    level = factor(level_labels[levs], levels = level_labels[levs]),
    n = as.integer(table(droplevels(dat[[variable]]))[levs]),
    fit = exp(p$fit),
    low = exp(p$fit - qnorm(0.975) * p$se.fit),
    high = exp(p$fit + qnorm(0.975) * p$se.fit)
  )
  ggplot(est, aes(x = level, y = fit)) +
    geom_pointrange(aes(ymin = low, ymax = high), colour = col_main, size = 0.3, linewidth = 0.6) +
    geom_text(aes(y = high, label = paste0("n = ", n)), vjust = -0.6, size = 2.3, colour = "grey30") +
    scale_x_discrete(labels = \(x) str_replace_all(x, c("/" = "/\n", " \\+ " = " +\n"))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = x_label, y = text_labels[[lang]]$density)
}

time_of_day_panel <- function(id, lang = "en") {
  lab <- text_labels[[lang]]
  factor_panel(id, "time_of_day_2", lab$time_of_day_levels, lab$time_of_day, lang)
}

loudness_class_panel <- function(id, lang = "en") {
  factor_panel(id, "loudness_class", text_labels[[lang]]$loudness_levels,
               axis_labels[[lang]][axis_labels$term == "loudness_s"], lang)
}

tag_panels <- function(p, title = NULL) {
  p + plot_annotation(tag_levels = "A", title = title,
                      theme = theme(plot.title = element_text(size = 9, face = "bold")))
}

#----------------------------------------------------------#
# Figures 1-5: hypothesis effects, English and Czech -----
#----------------------------------------------------------#
file_suffix <- c(en = "", cs = "_cs")

for (lang in c("en", "cs")) {
  lab <- text_labels[[lang]]
  save_lang <- \(plot, name, ...) save_figure(plot, paste0(name, file_suffix[[lang]]), ...)
  visit_1 <- labs(subtitle = lab$visit_1)
  visit_2 <- labs(subtitle = lab$visit_2)

  fig_1 <- tag_panels(
    effect_panel("H1", "spruce_s", lang) +
      effect_panel("H2", "ndop_s", lang, x_jitter = 0.15) +
      effect_panel("H3", "prob_s", lang)
  )
  save_lang(fig_1, "fig1_h1_h3_habitat", height_mm = 65)

  fig_2 <- tag_panels(
    effect_panel("H4", "calamity_s", lang) +
      effect_panel("H5", "dry_s", lang) +
      effect_panel("H6", "harvest_s", lang)
  )
  save_lang(fig_2, "fig2_h4_h6_disturbance", height_mm = 70)

  fig_3 <- tag_panels(
    effect_panel("H7", "doy_1_s", lang, x_jitter = 0.3) +
      effect_panel("H8", "doy_2_s", lang, x_jitter = 0.3) +
      time_of_day_panel("H8", lang)
  )
  save_lang(fig_3, "fig3_h7_h8_survey_date", height_mm = 70)

  # Spontaneous and playback-provoked detections separately
  fig_3b <- tag_panels(
    effect_panel("H7_spont", "doy_1_s", lang, x_jitter = 0.3) +
      effect_panel("H8_spont", "doy_2_s", lang, x_jitter = 0.3),
    title = lab$spontaneous
  )
  save_lang(fig_3b, "fig3b_h7_h8_survey_date_spontaneous", width_mm = width_double * 2 / 3, height_mm = 75)

  fig_3c <- tag_panels(
    effect_panel("H7_prov", "doy_1_s", lang, x_jitter = 0.3) +
      effect_panel("H8_prov", "doy_2_s", lang, x_jitter = 0.3),
    title = lab$provoked
  )
  save_lang(fig_3c, "fig3c_h7_h8_survey_date_provoked", width_mm = width_double * 2 / 3, height_mm = 75)

  # Without the 119-dB megaphone: values above 100 dB are uncertain, and four
  # fields at 119 dB stretched the axis and the fitted line
  fig_4 <- tag_panels(
    (effect_panel("H9_v1_no_mega", "loudness_s", lang, x_jitter = 0.3) + visit_1) +
      (effect_panel("H9_v2_no_mega", "loudness_s", lang, x_jitter = 0.3) + visit_2)
  )
  save_lang(fig_4, "fig4_h9_device_loudness", height_mm = 70)

  # Loudness classes, with everything >= 100 dB pooled
  fig_4b <- tag_panels(
    (loudness_class_panel("H9_v1_class", lang) + visit_1) +
      (loudness_class_panel("H9_v2_class", lang) + visit_2)
  )
  save_lang(fig_4b, "fig4b_h9_device_loudness_classes", height_mm = 70)

  fig_5 <- tag_panels(
    (effect_panel("H10_H11_pres", "lat_s", lang) + effect_panel("H10_H11_pres", "lon_s", lang)) /
      (effect_panel("H10_H11_dens", "lat_s", lang) + effect_panel("H10_H11_dens", "lon_s", lang))
  )
  save_lang(fig_5, "fig5_h10_h11_geography", height_mm = 120)

  fig_5b <- tag_panels(
    effect_panel("H10_H11_dens_selonly", "lat_s", lang) +
      effect_panel("H10_H11_dens_selonly", "lon_s", lang),
    title = lab$selected
  )
  save_lang(fig_5b, "fig5b_h10_h11_geography_selected", width_mm = width_double * 2 / 3, height_mm = 70)

  fig_5c <- tag_panels(
    effect_panel("H10_H11_dens_selonly_lowcal", "lat_s", lang) +
      effect_panel("H10_H11_dens_selonly_lowcal", "lon_s", lang),
    title = lab$selected_low_calamity
  )
  save_lang(fig_5c, "fig5c_h10_h11_geography_selected_low_calamity", width_mm = width_double * 2 / 3, height_mm = 70)
}

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
  scale_fill_gradient(
    name = density_label, low = pal_sequential[1], high = pal_sequential[2],
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
    # Brown = fewer, green = more detections than predicted
    scale_fill_distiller(
      name = "Scaled residual (DHARMa)", palette = "BrBG", direction = 1,
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

#----------------------------------------------------------#
# Supplementary: effects with and without selection status -----
#----------------------------------------------------------#
selection_pairs <- tribble(
  ~main,          ~adjusted,          ~term,        ~name,
  "H1",           "H1_sel",           "spruce_s",   "H1 – Spruce forest area",
  "H2",           "H2_sel",           "ndop_s",     "H2 – NDOP records",
  "H3",           "H3_sel",           "prob_s",     "H3 – Modelled occurrence probability",
  "H4",           "H4_sel",           "calamity_s", "H4 – Bark-beetle calamity",
  "H5",           "H5_sel",           "dry_s",      "H5 – Standing dead spruce",
  "H6",           "H6_sel",           "harvest_s",  "H6 – Salvage logging",
  "H10_H11_dens", "H10_H11_dens_sel", "lat_s",      "H10 – Latitude",
  "H10_H11_dens", "H10_H11_dens_sel", "lon_s",      "H11 – Longitude"
)

final_estimate <- function(id, term) {
  s <- spatial %>% filter(model_id == id)
  est <- if (s$use_spatial) s$estimates_after[[1]] else s$estimates_before[[1]]
  est %>% filter(term == !!term) %>% select(-term)
}

selection_effects <- selection_pairs %>%
  pivot_longer(c(main, adjusted), names_to = "fit", values_to = "model_id") %>%
  mutate(est = map2(model_id, term, final_estimate)) %>%
  unnest(est) %>%
  left_join(predictor_info %>% select(term, step_label), by = "term") %>%
  mutate(
    label = paste0(name, "\n(", step_label, ")"),
    label = forcats::fct_rev(forcats::fct_inorder(label)),
    fit = factor(fit, levels = c("main", "adjusted"),
                 labels = c("Main model", "Adjusted for selection status"))
  )

fig_s2 <- ggplot(selection_effects, aes(x = exp(estimate), y = label, colour = fit, shape = fit)) +
  geom_vline(xintercept = 1, colour = "grey50", linewidth = 0.3) +
  geom_pointrange(aes(xmin = exp(conf.low), xmax = exp(conf.high)),
                  position = position_dodge(width = 0.55), size = 0.3, linewidth = 0.5) +
  scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1, 2)) +
  scale_colour_manual(values = c(col_main, col_accent), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  labs(x = "Incidence rate ratio (95% CI, log scale)", y = NULL) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25))
save_figure(fig_s2, "figS2_selection_sensitivity", width_mm = width_single * 1.4, height_mm = 110)

message("04: figures written to ", dir_figures)
