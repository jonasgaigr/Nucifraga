#----------------------------------------------------------#
# 00 Setup: packages, paths, constants, shared plot theme -----
#----------------------------------------------------------#
# Sourced at the top of every numbered script so each one can also run on its
# own. Packages are not installed here; run `renv::restore()` once instead so
# that everyone uses the versions pinned in renv.lock.

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(spdep)
  library(glmmTMB)
  library(DHARMa)
  library(broom.mixed)
})

#----------------------------------------------------------#
# Paths -----
#----------------------------------------------------------#
path_raw <- here("Data", "Input", "Oresnik_statistika.xlsx")
dir_processed <- here("Data", "Processed")
dir_models <- here("Outputs", "Models")
dir_figures <- here("Outputs", "Figures")
dir_tables <- here("Outputs", "Tables")

walk(c(dir_processed, dir_models, dir_figures, dir_tables),
     dir.create, recursive = TRUE, showWarnings = FALSE)

path_clean_rds <- file.path(dir_processed, "fields_clean.rds")
path_models_rds <- file.path(dir_models, "models.rds")
path_spatial_rds <- file.path(dir_models, "spatial.rds")

#----------------------------------------------------------#
# Constants -----
#----------------------------------------------------------#
seed <- 20251010
set.seed(seed)

# S-JTSK / Krovak East-North. Verified in 01_import_clean.R against the KFME
# grid code of every field, so a silent axis swap would stop the pipeline.
crs_sjtsk <- 5514
crs_wgs84 <- 4326

# Protocol window for transect length (km)
line_protocol_km <- c(4, 6)

# Transect length at which detection probability is shown in figures; the
# midpoint of the protocol window keeps the presence panels interpretable.
line_reference_km <- 5

# Survey windows from the protocol; visits outside are treated as entry errors
# (e.g. visit 1 dated after visit 2) and excluded from the date models.
window_visit_1 <- 2:4 # February-April
window_visit_2 <- 6:7 # June-July

# Unit of the counts: territories (confirmed by the data owner; the source file
# itself does not say). Every axis label and reporting sentence follows this.
count_unit <- "territories"

# Predictors are centred and scaled so that the reported IRR/OR refers to a
# biologically meaningful step and the spatial models converge reliably.
# Only the slope interpretation matters; centring does not change it.
predictor_info <- tribble(
  ~term,        ~raw,            ~centre, ~scale, ~step_label,                 ~axis_label,                              ~phrase,
  "spruce_s",   "spruce_km2",    0,       1,      "per 1 km² of spruce forest", "Spruce forest area (km²)",               "spruce forest area",
  "ndop_s",     "ndop_records",  0,       1,      "per 1 NDOP record",          "NDOP records (n)",                       "the number of NDOP records",
  "prob_s",     "prob_occ",      0,       0.1,    "per 0.1 increase",           "Modelled occurrence probability",        "modelled occurrence probability",
  "calamity_s", "prop_calamity", 0,       0.1,    "per 10 percentage points",   "Bark-beetle calamity area (proportion)", "the proportion of bark-beetle calamity area",
  "dry_s",      "prop_dry",      0,       0.1,    "per 10 percentage points",   "Standing dead spruce area (proportion)", "the proportion of standing dead spruce",
  "harvest_s",  "prop_harvest",  0,       0.1,    "per 10 percentage points",   "Salvage-logged area (proportion)",       "the proportion of salvage-logged area",
  "doy_1_s",    "doy_1",         80,      10,     "per 10 days",                "Date of visit 1",                        "the date of visit 1",
  "doy_2_s",    "doy_2",         190,     10,     "per 10 days",                "Date of visit 2",                        "the date of visit 2",
  "loudness_s", "loudness_db",   85,      10,     "per 10 dB",                  "Device loudness (dB at 1 m)",            "device loudness",
  "lat_s",      "lat",           50,      1,      "per 1° latitude (~111 km)",  "Latitude (°N)",                          "latitude",
  "lon_s",      "lon",           15,      1,      "per 1° longitude (~72 km)",  "Longitude (°E)",                         "longitude"
)

#----------------------------------------------------------#
# Shared model helpers -----
#----------------------------------------------------------#
converged <- function(fit) {
  fit$fit$convergence == 0 && isTRUE(fit$sdr$pdHess)
}

build_formula <- function(response, predictors, random = NA, extra = NULL) {
  rhs <- c(predictors, "offset(log(line_km))", if (!is.na(random)) random, extra)
  as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
}

# McFadden pseudo-R² of the fixed effects: the null model keeps the offset and
# any random/spatial terms, so only the predictors' contribution is measured.
pseudo_r2 <- function(fit, data, response, random = NA, extra = NULL) {
  null_fit <- update(fit, formula = build_formula(response, "1", random, extra), data = data)
  as.numeric(1 - logLik(fit) / logLik(null_fit))
}

# Exponential spatial covariance on the field centroids (km, EPSG:5514)
spatial_term <- "exp(pos + 0 | group)"

add_spatial_coords <- function(data) {
  data %>% mutate(pos = numFactor(x_km, y_km), group = factor(1))
}

# Refit with the spatial covariance, keeping family, zero-inflation and random
# effects. The spatial random effect also absorbs extra-Poisson variation, and
# in these data NB2 theta then runs to infinity, which leaves the NB2 fit
# unconverged or on the boundary. A Poisson model with the spatial random
# effect (a spatial Poisson-lognormal model) is then the same model, estimated
# properly.
fit_spatial <- function(fit, data, response, predictors, random = NA) {
  sp_data <- add_spatial_coords(data)
  sp_formula <- build_formula(response, predictors, random, spatial_term)
  fit_sp <- update(fit, formula = sp_formula, data = sp_data)
  if (family(fit_sp)$family == "nbinom2" && (!converged(fit_sp) || sigma(fit_sp) > 1e3)) {
    fit_sp <- update(fit, formula = sp_formula, data = sp_data, family = poisson)
  }
  fit_sp
}

#----------------------------------------------------------#
# Plot theme and palettes -----
#----------------------------------------------------------#
# Journal column widths (mm)
width_single <- 85
width_double <- 175

# Okabe-Ito: designed to stay distinguishable under the common colour-vision
# deficiencies.
col_main <- "#0072B2"
col_accent <- "#D55E00"
col_points <- "grey35"
pal_categorical <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")

theme_nuc <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
      axis.line = element_line(colour = "grey30", linewidth = 0.3),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.3),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1, colour = "grey20"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      plot.tag = element_text(size = base_size + 2, face = "bold"),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = margin(4, 6, 4, 4)
    )
}
theme_set(theme_nuc())

# Saves one figure as vector PDF plus 600-dpi PNG and TIFF.
save_figure <- function(plot, name, width_mm = width_double, height_mm = 70) {
  base <- file.path(dir_figures, name)
  ggsave(paste0(base, ".pdf"), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf)
  ggsave(paste0(base, ".png"), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(paste0(base, ".tiff"), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = 600, device = ragg::agg_tiff, bg = "white",
         compression = "lzw")
  invisible(base)
}
