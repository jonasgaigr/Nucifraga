#----------------------------------------------------------#
# 03 Spatial autocorrelation checks and spatial refits -----
#----------------------------------------------------------#
# Input:  Outputs/Models/models.rds
# Output: Outputs/Models/spatial.rds
#
# Every model's residuals are tested with DHARMa (inverse-distance weights),
# Moran's I on a symmetric 8-nearest-neighbour graph, and a distance-class
# correlogram. Where any of the three is significant (for the correlogram: the
# first, 0-10 km class) the model is refitted with an
# exponential spatial covariance on the field centroids (glmmTMB). The fixed
# effects stay as they are: for H10/H11 latitude and longitude are the
# hypotheses themselves, so only the residual dependence is modelled and the
# linear gradient is never absorbed into a smooth of the coordinates.

source(here::here("R", "00_setup.R"))

models <- readRDS(path_models_rds)

alpha_sac <- 0.05
k_neighbours <- 8
correlogram_breaks <- seq(0, 150, by = 10)
n_sim <- 1000

#----------------------------------------------------------#
# Helpers -----
#----------------------------------------------------------#
# Residuals conditional on the fitted random effects. DHARMa cannot simulate
# conditionally for glmmTMB spatial covariance structures, and unconditional
# residuals would keep the very spatial (or observer) pattern the random
# effect accounts for. Simulating from the conditional mean treats every model
# the same way; for fixed-effect-only models it equals DHARMa's default.
simulate_conditional <- function(fit) {
  set.seed(seed)
  fam <- family(fit)$family
  mu <- predict(fit, type = "conditional")
  p_zero <- if (length(fixef(fit)$zi) > 0) predict(fit, type = "zprob") else 0
  y_obs <- model.response(model.frame(fit))
  n <- length(mu)

  sims <- replicate(n_sim, {
    y <- switch(fam,
      nbinom2 = rnbinom(n, size = sigma(fit), mu = mu),
      poisson = rpois(n, mu),
      binomial = rbinom(n, 1, mu)
    )
    y[runif(n) < p_zero] <- 0
    y
  })

  createDHARMa(sims, observedResponse = y_obs, fittedPredictedResponse = mu,
               integerResponse = TRUE, seed = seed)
}

neighbour_weights <- function(coords) {
  knearneigh(as.matrix(coords), k = k_neighbours) %>%
    knn2nb() %>%
    make.sym.nb() %>%
    nb2listw(style = "W")
}

sac_tests <- function(sim, coords, corr) {
  dharma <- testSpatialAutocorrelation(sim, x = coords$x_km, y = coords$y_km,
                                       plot = FALSE)
  moran <- moran.test(residuals(sim), neighbour_weights(coords),
                      randomisation = TRUE)
  tibble(
    dharma_I = unname(dharma$statistic["observed"]),
    dharma_p = dharma$p.value,
    moran_I = unname(moran$estimate["Moran I statistic"]),
    moran_p = moran$p.value,
    moran_0_10_I = corr$moran_I[1],
    moran_0_10_p = corr$p_value[1]
  )
}

correlogram <- function(sim, coords) {
  xy <- as.matrix(coords)
  res <- residuals(sim)
  map_dfr(seq_len(length(correlogram_breaks) - 1), \(i) {
    nb <- dnearneigh(xy, correlogram_breaks[i], correlogram_breaks[i + 1])
    if (sum(card(nb) > 0) < 10) return(NULL)
    mt <- moran.test(res, nb2listw(nb, style = "W", zero.policy = TRUE),
                     zero.policy = TRUE, randomisation = TRUE)
    tibble(
      distance_km = mean(correlogram_breaks[i:(i + 1)]),
      moran_I = unname(mt$estimate["Moran I statistic"]),
      p_value = mt$p.value
    )
  })
}

# Same zero-inflation and random effects as the non-spatial fit. The spatial
# random effect also absorbs extra-Poisson variation, and in these data NB2
# theta then runs to infinity, which leaves the NB2 fit unconverged or on the
# boundary. A Poisson model with the spatial random effect (a spatial
# Poisson-lognormal model) is then the same model, estimated properly.
fit_spatial <- function(fit, data, response, predictors, random) {
  sp_data <- add_spatial_coords(data)
  sp_formula <- build_formula(response, predictors, random, spatial_term)
  fit_sp <- update(fit, formula = sp_formula, data = sp_data)
  if (family(fit_sp)$family == "nbinom2" && (!converged(fit_sp) || sigma(fit_sp) > 1e3)) {
    fit_sp <- update(fit, formula = sp_formula, data = sp_data, family = poisson)
  }
  fit_sp
}

spatial_family_label <- function(fit) {
  link <- family(fit)$link
  switch(family(fit)$family,
    nbinom2 = "negative binomial (NB2), log link",
    poisson = "Poisson, log link",
    binomial = paste0("binomial, ", link, " link")
  ) %>%
    paste("+ exponential spatial covariance")
}

# Range parameter of the exponential correlation exp(-d / range), in km.
# Read from the fitted correlation matrix, whose rows follow the factor levels
# of `pos` (not the data order), so coordinates come from the level names.
# Correlation drops to 0.05 at about 3 x range.
spatial_range_km <- function(fit) {
  cor_mat <- attr(VarCorr(fit)$cond$group, "correlation")
  d <- as.matrix(dist(parseNumLevels(rownames(cor_mat))))
  keep <- upper.tri(cor_mat) & cor_mat > 0.01 & cor_mat < 0.99
  median(-d[keep] / log(cor_mat[keep]))
}

tidy_fixed <- function(fit) {
  tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    filter(component == "cond", term != "(Intercept)")
}

#----------------------------------------------------------#
# Tests on non-spatial models -----
#----------------------------------------------------------#
spatial <- models %>%
  select(model_id, hypothesis, role, response, predictors, random, fit, model_data) %>%
  mutate(
    coords = map(model_data, \(d) select(d, x_km, y_km)),
    sim_before = map(fit, simulate_conditional),
    correlogram_before = map2(sim_before, coords, correlogram),
    sac_before = pmap(list(sim_before, coords, correlogram_before), sac_tests),
    # Global tests average over all neighbour pairs and can miss dependence
    # that is confined to adjacent fields, so the first distance class counts too
    sac_present = map_lgl(sac_before, \(t) t$dharma_p < alpha_sac | t$moran_p < alpha_sac |
                            t$moran_0_10_p < alpha_sac)
  )

#----------------------------------------------------------#
# Spatial refits where needed -----
#----------------------------------------------------------#
spatial <- spatial %>%
  mutate(
    fit_sp = pmap(list(fit, model_data, response, predictors, random, sac_present, model_id),
                  \(f, d, r, p, re, sac, id) {
                    if (!sac) return(NULL)
                    message("Spatial refit: ", id)
                    fit_spatial(f, d, r, p, re)
                  }),
    sp_converged = map_lgl(fit_sp, \(f) !is.null(f) && converged(f)),
    sim_after = map(fit_sp, \(f) if (is.null(f)) NULL else simulate_conditional(f)),
    correlogram_after = map2(sim_after, coords, \(s, c) if (is.null(s)) NULL else correlogram(s, c)),
    sac_after = pmap(list(sim_after, coords, correlogram_after),
                     \(s, c, k) if (is.null(s)) NULL else sac_tests(s, c, k)),
    sp_family_label = map_chr(fit_sp, \(f) if (is.null(f)) NA_character_ else spatial_family_label(f)),
    sp_range_km = map_dbl(fit_sp, \(f) if (is.null(f)) NA_real_ else spatial_range_km(f)),
    sp_sd = map_dbl(fit_sp, \(f) if (is.null(f)) NA_real_ else attr(VarCorr(f)$cond$group, "stddev")[1]),
    sp_aic = map_dbl(fit_sp, \(f) if (is.null(f)) NA_real_ else AIC(f)),
    sp_pseudo_r2 = pmap_dbl(list(fit_sp, model_data, response, random),
                            \(f, d, r, re) {
                              if (is.null(f)) return(NA_real_)
                              pseudo_r2(f, add_spatial_coords(d), r, re, spatial_term)
                            }),
    # Final model: the spatial fit wherever it was needed and converged
    use_spatial = sac_present & sp_converged,
    fit_final = pmap(list(use_spatial, fit_sp, fit), \(u, a, b) if (u) a else b),
    sim_final = pmap(list(use_spatial, sim_after, sim_before), \(u, a, b) if (u) a else b),
    estimates_before = map(fit, tidy_fixed),
    estimates_after = map(fit_sp, \(f) if (is.null(f)) NULL else tidy_fixed(f))
  )

saveRDS(spatial, path_spatial_rds)

spatial %>%
  transmute(
    model_id,
    sac_present,
    use_spatial,
    moran_before = map_dbl(sac_before, "moran_p"),
    dharma_before = map_dbl(sac_before, "dharma_p"),
    first_class_before = map_dbl(sac_before, "moran_0_10_p"),
    moran_after = map_dbl(sac_after, \(t) if (is.null(t)) NA_real_ else t$moran_p),
    dharma_after = map_dbl(sac_after, \(t) if (is.null(t)) NA_real_ else t$dharma_p),
    sp_range_km, sp_sd
  ) %>%
  print(n = Inf, width = Inf)

message("03: spatial checks done; ", sum(spatial$use_spatial), " of ",
        nrow(spatial), " models refitted with spatial covariance.")
