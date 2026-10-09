#----------------------------------------------------------#
# 05 Result tables and reporting text -----
#----------------------------------------------------------#
# Input:  Data/Processed/fields_clean.rds, Outputs/Models/{models,spatial}.rds
# Output: Outputs/Tables/model_results.{csv,xlsx}, model_diagnostics.csv,
#         descriptives.csv, Outputs/results_summary.md

source(here::here("R", "00_setup.R"))

fields <- readRDS(path_clean_rds)
models <- readRDS(path_models_rds)
spatial <- readRDS(path_spatial_rds)

fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
fmt_p <- function(p) if_else(p < 0.001, "< 0.001", paste("=", formatC(p, format = "f", digits = 3)))

term_labels <- predictor_info %>%
  select(term, term_label = axis_label, step_label, phrase) %>%
  bind_rows(tibble(
    term = c("time_of_day_2afternoon/evening", "time_of_day_2morning + evening"),
    term_label = c("Time of day: afternoon/evening (vs morning)", "Time of day: morning + evening (vs morning)"),
    step_label = "vs morning",
    phrase = c("afternoon/evening visits", "visits spanning morning and evening")
  ))

capitalise <- function(x) paste0(toupper(str_sub(x, 1, 1)), str_sub(x, 2))

term_hypothesis <- c(lat_s = "H10", lon_s = "H11")

#----------------------------------------------------------#
# Model results: one row per term, per model, per fit type -----
#----------------------------------------------------------#
sac_columns <- function(tests, suffix) {
  tests %>%
    map(\(t) if (is.null(t)) tibble(dharma_I = NA, dharma_p = NA, moran_I = NA, moran_p = NA,
                                    moran_0_10_I = NA, moran_0_10_p = NA) else t) %>%
    bind_rows() %>%
    rename_with(\(x) paste0("sac_", x, "_", suffix))
}

model_meta <- models %>%
  select(model_id, hypothesis, role, response_label, family_label, n, expected, focal,
         theta, aic, pseudo_r2) %>%
  left_join(spatial %>% select(model_id, sac_present, use_spatial, sp_converged, sp_family_label,
                               sp_range_km, sp_sd, sp_aic, sp_pseudo_r2),
            by = "model_id") %>%
  bind_cols(sac_columns(spatial$sac_before, "before"),
            sac_columns(spatial$sac_after, "after"))

estimates <- spatial %>%
  select(model_id, estimates_before, estimates_after) %>%
  pivot_longer(-model_id, names_to = "fit_type", values_to = "est") %>%
  filter(!map_lgl(est, is.null)) %>%
  unnest(est) %>%
  mutate(fit_type = recode(fit_type, estimates_before = "non-spatial", estimates_after = "spatial"))

classify_support <- function(expected, estimate, low, high) {
  clear <- low > 0 | high < 0
  sign <- if_else(estimate > 0, "+", "-")
  case_when(
    !clear ~ "not supported (95% CI includes no effect)",
    expected == "+/-" ~ paste0("supported (", if_else(sign == "+", "positive", "negative"), " effect)"),
    sign == expected ~ "supported",
    TRUE ~ "not supported (effect in the opposite direction)"
  )
}

model_results <- estimates %>%
  left_join(model_meta, by = "model_id") %>%
  left_join(term_labels, by = "term") %>%
  mutate(
    is_final = (fit_type == "spatial") == use_spatial,
    is_focal = map2_lgl(term, focal, \(t, f) t %in% str_split_1(f, ",")),
    term_hypothesis = if_else(is_focal, coalesce(term_hypothesis[term], hypothesis), "covariate"),
    expected_term = map2_chr(term, seq_along(term), \(t, i) {
      exps <- str_split_1(expected[i], ",")
      focals <- str_split_1(focal[i], ",")
      if (t %in% focals) exps[match(t, focals)] else NA_character_
    }),
    effect_type = if_else(str_detect(family_label, "cloglog"),
                          "rate ratio (cloglog)", "IRR"),
    support = if_else(is_focal & is_final,
                      classify_support(expected_term, estimate, conf.low, conf.high),
                      NA_character_),
    # Model-level statistics belong to the fit the row comes from
    aic_fit = if_else(fit_type == "spatial", sp_aic, aic),
    pseudo_r2_fit = if_else(fit_type == "spatial", sp_pseudo_r2, pseudo_r2)
  ) %>%
  transmute(
    model_id, hypothesis = term_hypothesis, role, fit_type, is_final,
    response = response_label,
    family_link = if_else(fit_type == "spatial", sp_family_label, family_label),
    n, term, term_label, phrase, step = step_label,
    estimate_log = estimate, se = std.error,
    ci_low_log = conf.low, ci_high_log = conf.high,
    effect_type,
    effect = exp(estimate), effect_ci_low = exp(conf.low), effect_ci_high = exp(conf.high),
    z = statistic, p_value = p.value,
    expected_sign = expected_term, support,
    theta_nb2 = if_else(fit_type == "spatial" & !str_detect(family_link, "^negative"), NA_real_, theta),
    aic = aic_fit, pseudo_r2_mcfadden = pseudo_r2_fit,
    spatial_range_km = if_else(fit_type == "spatial", sp_range_km, NA_real_),
    spatial_sd = if_else(fit_type == "spatial", sp_sd, NA_real_),
    sac_moran_I_before, sac_moran_p_before, sac_moran_0_10_I_before, sac_moran_0_10_p_before,
    sac_dharma_p_before,
    sac_moran_I_after, sac_moran_p_after, sac_moran_0_10_I_after, sac_moran_0_10_p_after,
    sac_dharma_p_after
  ) %>%
  arrange(factor(model_id, levels = models$model_id), desc(fit_type == "non-spatial"))

write_csv(model_results, file.path(dir_tables, "model_results.csv"))
writexl::write_xlsx(list(model_results = model_results), file.path(dir_tables, "model_results.xlsx"))

#----------------------------------------------------------#
# Diagnostics per model -----
#----------------------------------------------------------#
model_diagnostics <- models %>%
  select(model_id, n, family_label) %>%
  bind_cols(bind_rows(models$diagnostics)) %>%
  left_join(
    models %>%
      select(model_id, family_comparison) %>%
      unnest(family_comparison) %>%
      mutate(family = paste0("aic_", str_replace_all(tolower(family), "[^a-z0-9]+", "_"))) %>%
      pivot_wider(names_from = family, values_from = aic),
    by = "model_id"
  )
write_csv(model_diagnostics, file.path(dir_tables, "model_diagnostics.csv"))

#----------------------------------------------------------#
# Collinearity of latitude and longitude -----
#----------------------------------------------------------#
# VIF from the correlation of the coefficient estimates (as car::vif does);
# with two predictors this equals 1 / (1 - r^2).
vif_two <- function(fit) {
  v <- vcov(fit)$cond[c("lat_s", "lon_s"), c("lat_s", "lon_s")]
  r <- cov2cor(v)[1, 2]
  1 / (1 - r^2)
}

total_data <- models$model_data[[which(models$model_id == "H10_H11_dens")]]
geo_collinearity <- tibble(
  cor_x_y_sjtsk = cor(total_data$x_sjtsk, total_data$y_sjtsk),
  cor_lon_lat = cor(total_data$lon, total_data$lat),
  vif_presence = vif_two(models$fit[[which(models$model_id == "H10_H11_pres")]]),
  vif_density = vif_two(models$fit[[which(models$model_id == "H10_H11_dens")]])
)

#----------------------------------------------------------#
# Descriptives -----
#----------------------------------------------------------#
describe <- function(x, variable, group = "") {
  tibble(
    group = group, variable = variable,
    n = sum(!is.na(x)), n_missing = sum(is.na(x)),
    mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE), q25 = quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    median = median(x, na.rm = TRUE), q75 = quantile(x, 0.75, na.rm = TRUE, names = FALSE),
    max = max(x, na.rm = TRUE), n_zero = sum(x == 0, na.rm = TRUE),
    prop_zero = mean(x == 0, na.rm = TRUE)
  )
}

count_tbl <- function(variable, values, group) {
  tibble(group = group, variable = variable, n = as.integer(values))
}

descriptives <- bind_rows(
  count_tbl("fields in source file", nrow(fields), "sample size"),
  count_tbl("fields with transect length", sum(!is.na(fields$line_km)), "sample size"),
  count_tbl("fields with both visits", sum(fields$visit_2_done), "sample size"),
  count_tbl("fields in total-count models", nrow(total_data), "sample size"),
  models %>% transmute(group = "sample size", variable = paste("n in model", model_id), n),
  fields %>% count(line_class) %>% transmute(group = "effort: transect length class", variable = line_class, n),
  describe(fields$line_km, "transect length (km)", "effort"),
  describe(fields$count_total, "total count (max of visits)", "raw counts"),
  describe(fields$count_1, "count visit 1", "raw counts"),
  describe(fields$count_2, "count visit 2", "raw counts"),
  describe(fields$spont_1, "spontaneous count visit 1", "raw counts"),
  describe(fields$prov_1, "provoked count visit 1", "raw counts"),
  describe(fields$spont_2, "spontaneous count visit 2", "raw counts"),
  describe(fields$prov_2, "provoked count visit 2", "raw counts"),
  describe(fields$presence, "presence (0/1)", "raw counts"),
  describe(fields$density, "total count per km", "density"),
  describe(fields$count_1 / fields$line_km, "visit 1 count per km", "density"),
  describe(fields$count_2 / fields$line_km, "visit 2 count per km", "density"),
  describe(fields$loudness_db, "device loudness (dB at 1 m)", "covariates"),
  describe(fields$doy_1, "visit 1 day of year", "covariates"),
  describe(fields$doy_2, "visit 2 day of year", "covariates"),
  describe(fields$lat, "latitude (°N)", "coordinates"),
  describe(fields$lon, "longitude (°E)", "coordinates"),
  geo_collinearity %>%
    pivot_longer(everything(), names_to = "variable", values_to = "mean") %>%
    mutate(group = "collinearity (total-count data)")
)
write_csv(descriptives, file.path(dir_tables, "descriptives.csv"))

#----------------------------------------------------------#
# Reporting text -----
#----------------------------------------------------------#
final_rows <- model_results %>% filter(is_final, role == "main")

effect_sentence <- function(id, term_id, outcome_text) {
  r <- final_rows %>% filter(model_id == id, term == term_id)
  direction <- if (r$ci_low_log > 0) "increased" else if (r$ci_high_log < 0) "decreased" else "did not change clearly"
  model_text <- paste0(
    r$family_link, ", offset log(transect length), n = ", r$n, " fields"
  )
  glue_text <- paste0(
    outcome_text, " ", direction, " with ", r$phrase,
    " (", r$effect_type, " = ", fmt(r$effect), " ", r$step, ", 95% CI ",
    fmt(r$effect_ci_low), "–", fmt(r$effect_ci_high), ", p ", fmt_p(r$p_value),
    "; ", model_text, ")."
  )
  list(text = glue_text, support = r$support, row = r)
}

sac_text <- function(id) {
  s <- spatial %>% filter(model_id == id)
  before <- s$sac_before[[1]]
  base <- paste0("Residual spatial autocorrelation in the non-spatial fit: Moran's I = ",
                 fmt(before$moran_I, 3), " (p ", fmt_p(before$moran_p), "; 0–10 km class: I = ",
                 fmt(before$moran_0_10_I, 3), ", p ", fmt_p(before$moran_0_10_p), "), DHARMa p ",
                 fmt_p(before$dharma_p), ".")
  if (!s$sac_present) return(paste(base, "No spatial refit was needed."))
  if (!s$use_spatial) return(paste(base, "The spatial refit did not converge; the non-spatial estimate is reported."))
  after <- s$sac_after[[1]]
  ns <- s$estimates_before[[1]] %>% filter(term %in% str_split_1(models$focal[models$model_id == id], ","))
  sp <- s$estimates_after[[1]] %>% filter(term %in% ns$term)
  comparison <- paste0(
    "non-spatial ", ns$term, ": ", fmt(exp(ns$estimate)), " (", fmt(exp(ns$conf.low)), "–",
    fmt(exp(ns$conf.high)), ") vs spatial ", fmt(exp(sp$estimate)), " (",
    fmt(exp(sp$conf.low)), "–", fmt(exp(sp$conf.high)), ")",
    collapse = "; "
  )
  paste0(base, " The reported estimate comes from the refit with an exponential spatial covariance (range parameter ",
         fmt(s$sp_range_km, 1), " km, i.e. correlation < 0.05 beyond about ", fmt(3 * s$sp_range_km, 0),
         " km); residual Moran's I after the refit = ", fmt(after$moran_I, 3),
         " (p ", fmt_p(after$moran_p), "). Comparison of the focal effect: ", comparison, ".")
}

unit <- count_unit
density_text <- paste0("Nutcracker density (", unit, " per km of transect)")

summary_items <- list(
  list(h = "H1 – Spruce forest area", id = "H1", term = "spruce_s", out = density_text),
  list(h = "H2 – NDOP records", id = "H2", term = "ndop_s", out = density_text),
  list(h = "H3 – Modelled occurrence probability", id = "H3", term = "prob_s", out = density_text),
  list(h = "H4 – Bark-beetle calamity", id = "H4", term = "calamity_s", out = density_text),
  list(h = "H5 – Standing dead spruce", id = "H5", term = "dry_s", out = density_text),
  list(h = "H6 – Salvage logging", id = "H6", term = "harvest_s", out = density_text),
  list(h = "H7 – Date of visit 1", id = "H7", term = "doy_1_s", out = "Density detected in visit 1"),
  list(h = "H8 – Date of visit 2", id = "H8", term = "doy_2_s", out = "Density detected in visit 2"),
  list(h = "H9 – Device loudness, visit 1", id = "H9_v1", term = "loudness_s", out = "Density of provoked detections in visit 1"),
  list(h = "H9 – Device loudness, visit 2", id = "H9_v2", term = "loudness_s", out = "Density of provoked detections in visit 2"),
  list(h = "H10 – Latitude, presence (separate model)", id = "H10_pres", term = "lat_s", out = "The detection probability per field"),
  list(h = "H10 – Latitude, density (separate model)", id = "H10_dens", term = "lat_s", out = density_text),
  list(h = "H10 – Latitude, presence (joint with longitude)", id = "H10_H11_pres", term = "lat_s", out = "The detection probability per field"),
  list(h = "H10 – Latitude, density (joint with longitude)", id = "H10_H11_dens", term = "lat_s", out = density_text),
  list(h = "H11 – Longitude, presence (separate model)", id = "H11_pres", term = "lon_s", out = "The detection probability per field"),
  list(h = "H11 – Longitude, density (separate model)", id = "H11_dens", term = "lon_s", out = density_text),
  list(h = "H11 – Longitude, presence (joint with latitude)", id = "H10_H11_pres", term = "lon_s", out = "The detection probability per field"),
  list(h = "H11 – Longitude, density (joint with latitude)", id = "H10_H11_dens", term = "lon_s", out = density_text)
)

sensitivity_line <- function(id, label) {
  r <- model_results %>% filter(model_id == id, term == "loudness_s", is_final)
  paste0(label, ": ", r$effect_type, " = ", fmt(r$effect), " (95% CI ", fmt(r$effect_ci_low), "–",
         fmt(r$effect_ci_high), ", p ", fmt_p(r$p_value), ", n = ", r$n, ")")
}

h8_tod <- final_rows %>% filter(model_id == "H8", str_detect(term, "time_of_day"))
west_note <- "Longitude effects are per 1° towards the east; a ratio below 1 means more nutcrackers further west, which is what H11 predicts."

paragraphs <- map_chr(summary_items, \(it) {
  s <- effect_sentence(it$id, it$term, it$out)
  extra <- case_when(
    it$id == "H8" ~ paste0(
      " Time of day was included as a covariate (", paste0(h8_tod$term_label, ": IRR = ", fmt(h8_tod$effect),
      ", 95% CI ", fmt(h8_tod$effect_ci_low), "–", fmt(h8_tod$effect_ci_high), collapse = "; "), ")."),
    it$id == "H9_v1" ~ paste0(" Sensitivity: ", sensitivity_line("H9_v1_no_mega", "without the 119-dB megaphone"),
                              "; ", sensitivity_line("H9_v1_no_re", "without the observer random effect"), "."),
    it$id == "H9_v2" ~ paste0(" Time of day was included as a covariate. Sensitivity: ",
                              sensitivity_line("H9_v2_no_mega", "without the 119-dB megaphone"),
                              "; ", sensitivity_line("H9_v2_no_re", "without the observer random effect"), "."),
    str_detect(it$term, "lon") ~ paste0(" ", west_note),
    TRUE ~ ""
  )
  paste0("### ", it$h, "\n\n", s$text, extra, " **", capitalise(s$support), ".**\n\n",
         sac_text(it$id), "\n")
})

overview <- map_dfr(summary_items, \(it) {
  r <- final_rows %>% filter(model_id == it$id, term == it$term)
  tibble(
    Hypothesis = it$h,
    Effect = paste0(r$effect_type, " ", fmt(r$effect), " (", fmt(r$effect_ci_low), "–",
                    fmt(r$effect_ci_high), ") ", r$step),
    p = str_remove(fmt_p(r$p_value), "^= "),
    `Spatial model` = if_else(str_detect(r$family_link, "spatial"), "yes", "no"),
    Verdict = capitalise(r$support)
  )
})
overview_md <- c(
  paste0("| ", paste(names(overview), collapse = " | "), " |"),
  paste0("|", strrep("---|", ncol(overview))),
  pmap_chr(overview, \(...) paste0("| ", paste(c(...), collapse = " | "), " |"))
)

d_line <- descriptives %>% filter(variable == "transect length (km)")
line_counts <- fields %>% count(line_class) %>% tibble::deframe()

caveats <- c(
  paste0("**Count unit.** `Celkem_zjisteno_oresniku` equals the maximum of the two visit totals in every field; the source file does not say whether the counts are individuals or territories. Results are labelled '", unit, "' and must be relabelled if they are individuals."),
  paste0("**Effort.** Transect length ranged ", fmt(d_line$min, 1), "–", fmt(d_line$max, 1), " km (median ", fmt(d_line$median, 1), "); ", line_counts[["below 4 km"]], " fields were below and ", line_counts[["above 6 km"]], " above the 4–6 km protocol, ", line_counts[["missing"]], " had no length. The offset assumes detections scale proportionally with transect length; lengths > 10 km may be the sum of both visits (to verify)."),
  "**Total count is a maximum.** Because the total is the maximum of two visits, it is not strictly proportional to the line length per visit; fields with only one visit were excluded from total-count models.",
  "**Disturbance variables.** Calamity area is exactly dry + harvested area, so H4 is not independent of H5 and H6. Proportions use remaining spruce + calamity area as the pre-disturbance spruce extent (the old code divided by remaining spruce and dropped 8 fields with proportions > 1).",
  "**H9 is observational.** Loudness is a property of the device, and almost every observer used one device, so loudness is confounded with observer skill and with the region each observer covered. The observer random effect absorbs part of this; the 119-dB megaphone (4 fields) and 9 fields without loudness data are handled in sensitivity runs and exclusions.",
  "**Doba_2 is time of day, not provocation duration.** No duration variable exists in the data, so H9 visit 2 and H8 adjust for time of day instead.",
  paste0("**Geography as a proxy.** Latitude and longitude are weakly correlated (r = ", fmt(geo_collinearity$cor_lon_lat), ", VIF = ", fmt(geo_collinearity$vif_density), " in the joint density model), so the joint and separate models agree. The gradients are probably proxies for elevation, spruce forest cover and mountain ranges (border mountains in the north and west), not causal effects of position. No elevation or forest-cover data were added here."),
  "**Spatial confounding.** In the spatial refits the spatial random effect and the linear latitude/longitude terms compete for the same broad-scale variation. Attenuated or wider-CI gradients after the refit are expected and do not mean the non-spatial estimate was wrong; both are reported in `model_results.csv`.",
  "**Spatial count models are Poisson-lognormal.** In all spatial refits of the main count models the NB2 dispersion parameter ran to infinity, because the spatial random effect absorbs the extra-Poisson variation. Those refits therefore use a Poisson family with the spatial random effect, which is the boundary case of the same model (see `family_link`).",
  "**Conditional residuals.** Post-refit spatial tests use residuals conditional on the fitted spatial random effect (DHARMa cannot simulate conditionally for glmmTMB spatial structures). These remove autocorrelation partly by construction, so they confirm that the structure is adequate rather than providing an independent test.",
  "**Multiple testing.** Eleven hypotheses (18 focal tests) were evaluated without correction; individual p-values near 0.05 should be read with caution.",
  "**Sampling design.** `Vyber_2025` (selected vs not selected fields) differs strongly in mean count and was not modelled; if field selection depended on expected occurrence, the habitat effects may be partly design-driven."
)

legends <- c(
  "**Figure 1.** Effects of (A) spruce forest area, (B) number of NDOP records and (C) modelled occurrence probability on nutcracker density. Points: observed counts per km of transect (grey, semi-transparent; B jittered horizontally). Lines and ribbons: predictions and 95% CIs from the final models with a log(transect length) offset, shown per 1 km of transect. Square-root y axis.",
  "**Figure 2.** Effects of the proportion of (A) bark-beetle calamity, (B) standing dead spruce and (C) salvage-logged area on nutcracker density. Layout as in Figure 1.",
  "**Figure 3.** (A) Density detected in visit 1 vs date; (B) density detected in visit 2 vs date (time of day fixed at morning); (C) density in visit 2 by time of day at the mean date (n = number of fields).",
  "**Figure 4.** Density of provoked detections vs playback device loudness in (A) visit 1 and (B) visit 2. Predictions at population level (observer random effect = 0); visit 2 at morning.",
  paste0("**Figure 5.** Geographic gradients from the joint latitude + longitude models: (A, B) probability of detecting nutcrackers on a ", line_reference_km, "-km transect (cloglog model with offset), observed presence/absence jittered vertically; (C, D) density per km of transect."),
  "**Figure 6.** Nutcracker density per km of transect in each mapping field (KFME sixteenth cells; grey = transect length missing). Projection S-JTSK / Krovak East North (EPSG:5514); graticule in WGS84.",
  "**Figure 7.** DHARMa scaled residuals of the final (A) presence and (B) density models with latitude + longitude. Values near 0.5 are as expected; orange (near 0) = fewer and purple (near 1) = more detections than predicted.",
  "**Figure S1.** Moran's I correlograms (10-km distance classes) of the residuals of every main model before (blue) and, where needed, after (orange) adding an exponential spatial covariance. Filled points: p < 0.05."
)

summary_md <- c(
  "# Nutcracker analysis – results summary",
  "",
  paste0("Generated by `R/run_all.R` on ", Sys.Date(), ". All numbers come from `Outputs/Tables/model_results.csv` (rows with `is_final = TRUE`)."),
  "",
  "Effect sizes are incidence rate ratios (IRR) for density models (count with offset log(transect length), so the IRR acts on count per km). For presence models with a complementary log-log link and the same offset, exp(β) is a ratio of the underlying detection rates, interpretable like an IRR (not an odds ratio). Hypotheses count as supported when the 95% CI excludes 1 in the predicted direction.",
  "",
  "## Overview",
  "",
  overview_md,
  "",
  "## Hypotheses",
  "",
  paragraphs,
  "## Caveats",
  "",
  paste("-", caveats),
  "",
  "## Figure legends",
  "",
  paste0(legends, "\n")
)

write_lines(summary_md, here("Outputs", "results_summary.md"))

message("05: tables and results summary written.")
