#----------------------------------------------------------#
# 01 Import and clean the field data -----
#----------------------------------------------------------#
# Input:  Data/Input/Oresnik_statistika.xlsx (never modified)
# Output: Data/Processed/fields_clean.{rds,csv}, Outputs/Tables/data_checks.csv

source(here::here("R", "00_setup.R"))

#----------------------------------------------------------#
# Read -----
#----------------------------------------------------------#
# Everything as text first: several numeric columns contain "NO DATA" or "–",
# and dates are Excel serial numbers stored as text.
raw <- read_excel(path_raw, sheet = 1, col_types = "text")

# English snake_case names; mapping to the original columns is in the README.
fields <- raw %>%
  rename(
    field_id = ID,
    kfme_cell = POLE,
    spruce_m2 = SMRK_m2,
    ndop_records = Pocet_nalezu,
    prob_occ = prob_1h_mean,
    calamity_m2 = `Plocha kurovcove kalamity celkem (m2)`,
    dry_m2 = `Plocha_souse (m2)`,
    harvest_m2 = `Plocha_tezby (m2)`,
    selected_2025 = Vyber_2025,
    observer = Mapovatel,
    count_total = Celkem_zjisteno_oresniku,
    date_1 = Datum_1,
    spont_1 = Spontanne_1,
    prov_1 = Provokovano_1,
    count_1 = Termin_1,
    date_2 = Datum_2,
    time_of_day_2 = Doba_2,
    spont_2 = Spontanne_2,
    prov_2 = Provokovano_2,
    count_2 = Termin_2,
    device = Zarizeni,
    loudness_db = Hlasitost,
    line_raw = DELKA_LINII,
    y_sjtsk = Y,
    x_sjtsk = X
  )

missing_codes <- c("NO DATA", "–", "-", "")

excel_date <- function(x) {
  as_date(as.numeric(if_else(x %in% missing_codes, NA_character_, x)),
          origin = "1899-12-30")
}

as_num <- function(x) as.numeric(if_else(x %in% missing_codes, NA_character_, x))

fields <- fields %>%
  mutate(
    across(c(field_id, spruce_m2, ndop_records, prob_occ, calamity_m2, dry_m2,
             harvest_m2, count_total, spont_1, prov_1, count_1, spont_2,
             prov_2, count_2, loudness_db, line_raw, y_sjtsk, x_sjtsk),
           as_num),
    across(c(date_1, date_2), excel_date),
    across(c(field_id, ndop_records, count_total, spont_1, prov_1, count_1,
             spont_2, prov_2, count_2), as.integer)
  )

#----------------------------------------------------------#
# Structural checks: stop if the data no longer match the assumptions -----
#----------------------------------------------------------#
stopifnot(
  !anyDuplicated(fields$kfme_cell),
  all(fields$count_1 == fields$spont_1 + fields$prov_1),
  all(fields$count_2 == fields$spont_2 + fields$prov_2, na.rm = TRUE),
  # The total is the maximum over the two visits, not their sum
  all(fields$count_total == pmax(fields$count_1, fields$count_2))
)

#----------------------------------------------------------#
# Coordinates: verify CRS against the KFME grid code -----
#----------------------------------------------------------#
# X/Y are negative, so they are S-JTSK (EPSG:5514, East-North axis order):
# X = easting (more negative = further west), Y = northing (more negative =
# further south). The Krovak grid is rotated against true north, so the
# latitude/longitude hypotheses use WGS84 coordinates instead.
fields_sf <- st_as_sf(fields, coords = c("x_sjtsk", "y_sjtsk"),
                      crs = crs_sjtsk, remove = FALSE)
lonlat <- st_coordinates(st_transform(fields_sf, crs_wgs84))

# KFME cell centroid from its code: rows of 6' latitude from 56°N, columns of
# 10' longitude from 5°40'E; quadrant letters a = NW, b = NE, c = SW, d = SE.
kfme_centroid <- function(code) {
  row <- as.integer(str_sub(code, 1, 2))
  col <- as.integer(str_sub(code, 3, 4))
  q1 <- str_sub(code, 5, 5)
  q2 <- str_sub(code, 6, 6)
  south <- function(q) q %in% c("c", "d")
  east <- function(q) q %in% c("b", "d")
  tibble(
    lat_kfme = 56 - row / 10 - 0.05 * south(q1) - 0.025 * south(q2) - 0.0125,
    lon_kfme = 17 / 3 + col / 6 + east(q1) / 12 + east(q2) / 24 + 1 / 48
  )
}

fields <- fields %>%
  mutate(lon = lonlat[, "X"], lat = lonlat[, "Y"]) %>%
  bind_cols(kfme_centroid(.$kfme_cell)) %>%
  mutate(kfme_offset_m = sqrt(((lon - lon_kfme) * cos(lat * pi / 180) * 111320)^2 +
                                ((lat - lat_kfme) * 110570)^2))

if (max(fields$kfme_offset_m) > 50) {
  stop("Coordinates do not match KFME cell centroids; check CRS and axis order.")
}

#----------------------------------------------------------#
# Manual corrections documented in the source data -----
#----------------------------------------------------------#
corrections <- tribble(
  ~kfme_cell, ~variable, ~reason,
  "6048dc",   "prov_1",  "Device note: first visit done without playback, so no provoked count exists"
)

fields <- fields %>%
  mutate(
    prov_1 = if_else(kfme_cell %in% corrections$kfme_cell[corrections$variable == "prov_1"],
                     NA_integer_, prov_1),
    # One team is spelled in two ways
    observer = str_replace(observer, "Vítezslav", "Vítězslav")
  )

# Visit 1 entered in August, after visit 2, in two fields of one observer
# (confirmed entry errors). The month is read as 03 typed as 08: the corrected
# dates fall within this observer's other visit-1 dates (22 March-3 April) and
# on days on which the same observer surveyed neighbouring fields.
date_fixes <- tribble(
  ~kfme_cell, ~date_1_recorded, ~date_1_corrected,
  "5253ab",   "2025-08-30",     "2025-03-30",
  "5254bb",   "2025-08-23",     "2025-03-23"
) %>%
  mutate(across(starts_with("date"), as_date))

recorded <- fields$date_1[match(date_fixes$kfme_cell, fields$kfme_cell)]
if (!identical(recorded, date_fixes$date_1_recorded)) {
  stop("Source dates for the visit-1 fixes have changed; review date_fixes in 01_import_clean.R.")
}

fields <- fields %>%
  rows_update(date_fixes %>% select(kfme_cell, date_1 = date_1_corrected), by = "kfme_cell")

stopifnot(!any(fields$date_1 > fields$date_2, na.rm = TRUE))

#----------------------------------------------------------#
# Derived variables -----
#----------------------------------------------------------#
fields <- fields %>%
  mutate(
    # Visit 2 was not done where its components are missing, but count_2 was
    # entered as 0; a 0 would wrongly count as a survey with no detections.
    visit_2_done = !is.na(spont_2),
    count_2 = if_else(visit_2_done, count_2, NA_integer_),
    n_visits = 1L + visit_2_done,

    # Transect length is already in km (range 1.7-15.4). Guard against a
    # future export in metres.
    line_km = if (median(line_raw, na.rm = TRUE) > 100) line_raw / 1000 else line_raw,
    line_class = case_when(
      is.na(line_km) ~ "missing",
      line_km <= 0 ~ "zero",
      line_km < line_protocol_km[1] ~ "below 4 km",
      line_km > line_protocol_km[2] ~ "above 6 km",
      TRUE ~ "within 4-6 km"
    ),

    density = count_total / line_km,
    presence = as.integer(count_total > 0),

    spruce_km2 = spruce_m2 / 1e6,
    # SMRK_m2 is the spruce remaining after the calamity (confirmed by the data
    # owner). Calamity area equals dry + harvested area, and exceeds remaining
    # spruce in some fields. Remaining spruce + calamity is therefore the
    # pre-disturbance spruce extent, so the proportions stay within 0-1 and
    # no field has to be dropped.
    spruce_before_m2 = spruce_m2 + calamity_m2,
    prop_calamity = calamity_m2 / spruce_before_m2,
    prop_dry = dry_m2 / spruce_before_m2,
    prop_harvest = harvest_m2 / spruce_before_m2,

    doy_1 = yday(date_1),
    doy_2 = yday(date_2),
    in_window_1 = month(date_1) %in% window_visit_1,
    in_window_2 = month(date_2) %in% window_visit_2,

    # Doba_2 records time of day of visit 2, not provocation duration
    time_of_day_2 = factor(
      recode(time_of_day_2,
             "Ráno/dopoledne" = "morning",
             "Odpoledne/večer" = "afternoon/evening",
             "Ráno/dopoledne večer" = "morning + evening",
             "–" = NA_character_),
      levels = c("morning", "afternoon/evening", "morning + evening")
    ),

    # Fields selected for 2025 were most likely chosen by expected occurrence,
    # so selection status is a design covariate in sensitivity models
    selected_2025 = factor(recode(selected_2025, ano = "yes", ne = "no"), levels = c("no", "yes")),
    observer = factor(observer),
    x_km = x_sjtsk / 1000,
    y_km = y_sjtsk / 1000
  )

# Scaled model predictors, defined once in 00_setup.R
for (i in seq_len(nrow(predictor_info))) {
  p <- predictor_info[i, ]
  fields[[p$term]] <- (fields[[p$raw]] - p$centre) / p$scale
}

#----------------------------------------------------------#
# Flag records that need manual verification -----
#----------------------------------------------------------#
flag <- function(cond, issue) {
  fields %>%
    filter({{ cond }}) %>%
    transmute(kfme_cell, observer = as.character(observer), issue = issue)
}

data_checks <- bind_rows(
  flag(is.na(line_km), "Transect length missing ('NO DATA'); excluded from all models"),
  flag(line_km > 10, "Transect length > 10 km; confirmed as one visit's length (not both visits summed), kept"),
  flag(is.na(loudness_db), "Device loudness missing ('NO DATA'); excluded from H9"),
  flag(loudness_db >= 110, "Device loudness >= 110 dB (megaphone-type caller); H9 sensitivity run without it"),
  flag(str_detect(device, "musím ještě ověřit|nevím"), "Device identity marked uncertain by observer; verified by data owner, kept"),
  flag(kfme_cell %in% corrections$kfme_cell, "Visit 1 without playback; provoked count set to NA"),
  flag(!visit_2_done, "Visit 2 not done but Termin_2 coded 0; set to NA, field excluded from total-count models"),
  flag(is.na(date_1), "Date of visit 1 missing"),
  flag(is.na(date_2), "Date of visit 2 missing"),
  flag(!is.na(date_1) & !in_window_1, "Visit 1 outside February-April; excluded from H7"),
  flag(!is.na(date_2) & !in_window_2, "Visit 2 outside June-July; excluded from H8"),
  flag(kfme_cell %in% date_fixes$kfme_cell, "Visit 1 dated in August, after visit 2 (entry error); month corrected to March"),
  flag(spruce_m2 < calamity_m2, "Calamity area exceeds remaining spruce area; kept (proportions use remaining spruce + calamity)")
) %>%
  arrange(issue, kfme_cell)

#----------------------------------------------------------#
# Write -----
#----------------------------------------------------------#
fields_out <- fields %>%
  select(-line_raw, -lat_kfme, -lon_kfme, -kfme_offset_m)

saveRDS(fields_out, path_clean_rds)
write_csv(fields_out, file.path(dir_processed, "fields_clean.csv"))
write_csv(data_checks, file.path(dir_tables, "data_checks.csv"))

message(
  "01: ", nrow(fields_out), " fields cleaned; CRS EPSG:", crs_sjtsk,
  " verified (max offset from KFME centroid ", round(max(fields$kfme_offset_m), 2), " m); ",
  nrow(data_checks), " records flagged."
)
