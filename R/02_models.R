#----------------------------------------------------------#
# 02 Fit non-spatial models for H1-H11 -----
#----------------------------------------------------------#
# Input:  Data/Processed/fields_clean.rds
# Output: Outputs/Models/models.rds (fits + DHARMa diagnostics)
#
# All count models use the raw count with offset(log(line_km)), so
# coefficients act on density (count per km of transect) while the response
# keeps its count distribution. Fitted in glmmTMB throughout, so that Poisson,
# NB2, zero-inflated, random-effect and spatial variants share one engine and
# comparable AIC; nbinom2 is the same likelihood as MASS::glm.nb.

source(here::here("R", "00_setup.R"))

fields <- readRDS(path_clean_rds)

#----------------------------------------------------------#
# Analysis data sets -----
#----------------------------------------------------------#
# Fields with only one visit are excluded from total-count models: the total
# is the maximum over visits, so one visit is lower effort that the
# line-length offset cannot correct.
has_line <- fields %>% filter(!is.na(line_km))

data_sets <- list(
  total = has_line %>% filter(visit_2_done),
  visit_1 = has_line %>% filter(in_window_1),
  visit_2 = has_line %>% filter(visit_2_done, in_window_2, !is.na(time_of_day_2)),
  prov_1 = has_line %>% filter(!is.na(prov_1), !is.na(loudness_db)),
  prov_2 = has_line %>% filter(visit_2_done, !is.na(loudness_db), !is.na(time_of_day_2))
)
data_sets$prov_1_no_megaphone <- data_sets$prov_1 %>% filter(loudness_db < 110)
data_sets$prov_2_no_megaphone <- data_sets$prov_2 %>% filter(loudness_db < 110)

#----------------------------------------------------------#
# Model specifications -----
#----------------------------------------------------------#
# expected: hypothesised sign of the focal slope(s) on the log scale.
# For H11 ("further west"), west is lower longitude, so the sign is negative.
model_specs <- tribble(
  ~model_id,        ~hypothesis, ~type,      ~response,     ~predictors,                 ~focal,             ~expected, ~random,           ~data,                 ~role,
  "H1",             "H1",        "count",    "count_total", "spruce_s",                  "spruce_s",         "+",       NA,                "total",               "main",
  "H2",             "H2",        "count",    "count_total", "ndop_s",                    "ndop_s",           "+",       NA,                "total",               "main",
  "H3",             "H3",        "count",    "count_total", "prob_s",                    "prob_s",           "+",       NA,                "total",               "main",
  "H4",             "H4",        "count",    "count_total", "calamity_s",                "calamity_s",       "-",       NA,                "total",               "main",
  "H5",             "H5",        "count",    "count_total", "dry_s",                     "dry_s",            "-",       NA,                "total",               "main",
  "H6",             "H6",        "count",    "count_total", "harvest_s",                 "harvest_s",        "-",       NA,                "total",               "main",
  "H7",             "H7",        "count",    "count_1",     "doy_1_s",                   "doy_1_s",          "+/-",     NA,                "visit_1",             "main",
  "H8",             "H8",        "count",    "count_2",     "doy_2_s + time_of_day_2",   "doy_2_s",          "+/-",     NA,                "visit_2",             "main",
  "H9_v1",          "H9",        "count",    "prov_1",      "loudness_s",                "loudness_s",       "+",       "(1 | observer)",  "prov_1",              "main",
  "H9_v2",          "H9",        "count",    "prov_2",      "loudness_s + time_of_day_2","loudness_s",       "+",       "(1 | observer)",  "prov_2",              "main",
  "H10_pres",       "H10",       "presence", "presence",    "lat_s",                     "lat_s",            "+",       NA,                "total",               "main",
  "H11_pres",       "H11",       "presence", "presence",    "lon_s",                     "lon_s",            "-",       NA,                "total",               "main",
  "H10_H11_pres",   "H10/H11",   "presence", "presence",    "lat_s + lon_s",             "lat_s,lon_s",      "+,-",     NA,                "total",               "main",
  "H10_dens",       "H10",       "count",    "count_total", "lat_s",                     "lat_s",            "+",       NA,                "total",               "main",
  "H11_dens",       "H11",       "count",    "count_total", "lon_s",                     "lon_s",            "-",       NA,                "total",               "main",
  "H10_H11_dens",   "H10/H11",   "count",    "count_total", "lat_s + lon_s",             "lat_s,lon_s",      "+,-",     NA,                "total",               "main",
  # Sensitivity: without the 119-dB megaphone, and without the observer term
  "H9_v1_no_mega",  "H9",        "count",    "prov_1",      "loudness_s",                "loudness_s",       "+",       "(1 | observer)",  "prov_1_no_megaphone", "sensitivity",
  "H9_v2_no_mega",  "H9",        "count",    "prov_2",      "loudness_s + time_of_day_2","loudness_s",       "+",       "(1 | observer)",  "prov_2_no_megaphone", "sensitivity",
  "H9_v1_no_re",    "H9",        "count",    "prov_1",      "loudness_s",                "loudness_s",       "+",       NA,                "prov_1",              "sensitivity",
  "H9_v2_no_re",    "H9",        "count",    "prov_2",      "loudness_s + time_of_day_2","loudness_s",       "+",       NA,                "prov_2",              "sensitivity"
)

response_label <- c(
  count_total = "Total count (max of 2 visits)",
  count_1 = "Count, visit 1",
  count_2 = "Count, visit 2",
  prov_1 = "Provoked count, visit 1",
  prov_2 = "Provoked count, visit 2",
  presence = "Presence (detected in either visit)"
)

#----------------------------------------------------------#
# Fitting helpers -----
#----------------------------------------------------------#
dharma_summary <- function(sim) {
  tibble(
    dispersion_ratio = unname(testDispersion(sim, plot = FALSE)$statistic),
    dispersion_p = testDispersion(sim, plot = FALSE)$p.value,
    zi_ratio = unname(testZeroInflation(sim, plot = FALSE)$statistic),
    zi_p = testZeroInflation(sim, plot = FALSE)$p.value,
    uniformity_p = testUniformity(sim, plot = FALSE)$p.value,
    outlier_p = testOutliers(sim, type = "bootstrap", plot = FALSE)$p.value
  )
}

fit_spec <- function(model_id, type, response, predictors, random, data, ...) {
  dat <- data_sets[[data]]
  f <- build_formula(response, predictors, random)

  if (type == "presence") {
    final <- glmmTMB(f, family = binomial(link = "cloglog"), data = dat)
    sim <- simulateResiduals(final, n = 1000, seed = seed)
    return(list(
      fit = final, sim = sim, data = dat, formula = f,
      family_label = "binomial, cloglog link",
      family_comparison = tibble(family = "binomial (cloglog)", aic = AIC(final)),
      diagnostics = dharma_summary(sim),
      pseudo_r2 = pseudo_r2(final, dat, response, random)
    ))
  }

  fit_po <- glmmTMB(f, family = poisson, data = dat)
  fit_nb <- glmmTMB(f, family = nbinom2, data = dat)
  sim_po <- simulateResiduals(fit_po, n = 1000, seed = seed)
  sim_nb <- simulateResiduals(fit_nb, n = 1000, seed = seed)

  final <- fit_nb
  sim <- sim_nb
  family_label <- "negative binomial (NB2), log link"
  comparison <- tibble(family = c("Poisson", "NB2"), aic = c(AIC(fit_po), AIC(fit_nb)))

  # Once a random effect absorbs the extra variation, NB2 theta can run to
  # infinity (NB2 -> Poisson) and the fit fails or sits on the boundary; the
  # Poisson model is then the same model, estimated properly.
  if (!converged(fit_nb) || sigma(fit_nb) > 1e3) {
    final <- fit_po
    sim <- sim_po
    family_label <- "Poisson, log link (NB2 theta at boundary)"
  }

  # Zero-inflation only when the NB residuals call for it, and only kept if it
  # clearly improves the fit.
  zi_test <- testZeroInflation(sim, plot = FALSE)
  if (zi_test$p.value < 0.05) {
    fit_zinb <- update(fit_nb, ziformula = ~1, data = dat)
    comparison <- bind_rows(comparison, tibble(family = "ZINB", aic = AIC(fit_zinb)))
    if (AIC(fit_zinb) < AIC(final) - 2) {
      final <- fit_zinb
      sim <- simulateResiduals(final, n = 1000, seed = seed)
      family_label <- "zero-inflated NB2, log link"
    }
  }

  list(
    fit = final, sim = sim, data = dat, formula = f,
    family_label = family_label,
    family_comparison = comparison,
    diagnostics = dharma_summary(sim) %>%
      mutate(
        poisson_dispersion_ratio = unname(testDispersion(sim_po, plot = FALSE)$statistic),
        poisson_dispersion_p = testDispersion(sim_po, plot = FALSE)$p.value,
        poisson_zi_p = testZeroInflation(sim_po, plot = FALSE)$p.value
      ),
    pseudo_r2 = pseudo_r2(final, dat, response, random)
  )
}

#----------------------------------------------------------#
# Fit -----
#----------------------------------------------------------#
models <- model_specs %>%
  mutate(
    response_label = unname(response_label[response]),
    result = pmap(model_specs, fit_spec),
    fit = map(result, "fit"),
    sim = map(result, "sim"),
    model_data = map(result, "data"),
    formula = map(result, "formula"),
    family_label = map_chr(result, "family_label"),
    family_comparison = map(result, "family_comparison"),
    diagnostics = map(result, "diagnostics"),
    pseudo_r2 = map_dbl(result, "pseudo_r2"),
    n = map_int(model_data, nrow),
    aic = map_dbl(fit, AIC),
    # NB2 overdispersion parameter (theta = size); larger = closer to Poisson
    theta = map2_dbl(fit, family_label, \(f, lab) if (str_detect(lab, "^(negative|zero)")) sigma(f) else NA_real_)
  ) %>%
  select(-result)

saveRDS(models, path_models_rds)

models %>%
  select(model_id, n, family_label, aic, theta) %>%
  bind_cols(bind_rows(models$diagnostics)) %>%
  print(n = Inf, width = Inf)

message("02: ", nrow(models), " models fitted and saved to ", path_models_rds)
