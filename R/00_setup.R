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

# Unit of the counts is not stated in the source file (see README). Change it
# here and every axis label and reporting sentence follows.
count_unit <- "territories"

# Predictors are centred and scaled so that the reported IRR/OR refers to a
# biologically meaningful step and the spatial models converge reliably.
# Only the slope interpretation matters; centring does not change it.
predictor_info <- tribble(
  ~term,        ~raw,             ~centre, ~scale, ~step_label,                     ~axis_label,
  "spruce_s",   "spruce_km2",     0,       1,      "per 1 km² of spruce forest",    "Spruce forest area (km²)",
  "ndop_s",     "ndop_records",   0,       1,      "per 1 NDOP record",             "NDOP records (n)",
  "prob_s",     "prob_occ",       0,       0.1,    "per 0.1 increase",              "Modelled occurrence probability",
  "calamity_s", "prop_calamity",  0,       0.1,    "per 10 percentage points",      "Bark-beetle calamity area (proportion)",
  "dry_s",      "prop_dry",       0,       0.1,    "per 10 percentage points",      "Standing dead spruce area (proportion)",
  "harvest_s",  "prop_harvest",   0,       0.1,    "per 10 percentage points",      "Salvage-logged area (proportion)",
  "doy_1_s",    "doy_1",          80,      10,     "per 10 days",                   "Date of visit 1",
  "doy_2_s",    "doy_2",          190,     10,     "per 10 days",                   "Date of visit 2",
  "loudness_s", "loudness_db",    85,      10,     "per 10 dB",                     "Device loudness (dB at 1 m)",
  "lat_s",      "lat",            50,      1,      "per 1° latitude (~111 km)",     "Latitude (°N)",
  "lon_s",      "lon",            15,      1,      "per 1° longitude (~72 km)",     "Longitude (°E)"
)

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
