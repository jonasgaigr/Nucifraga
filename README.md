# Nutcracker (*Nucifraga caryocatactes*) survey analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Analysis of the 2025 Czech nutcracker mapping in KFME sixteenth-cells (about
3 × 2.8 km). Two visits per field record spontaneous and playback-provoked
detections along transects. This repository tests how nutcracker density
depends on habitat, disturbance, survey timing, playback device loudness and
geographic position. Survey effort (transect length) is corrected for
throughout, and residual spatial autocorrelation is checked for every model.

---

## Repository layout

```
Data/
├── Input/Oresnik_statistika.xlsx   # raw data – never modified
└── Processed/fields_clean.csv      # cleaned data written by R/01
R/
├── 00_setup.R          # packages, paths, constants, shared helpers and ggplot theme
├── 01_import_clean.R   # import, type conversion, CRS check, derived variables, data flags
├── 02_models.R         # Poisson vs NB2 (vs ZINB) models with effort offset, DHARMa diagnostics
├── 03_spatial.R        # spatial autocorrelation tests, spatial refits where needed
├── 04_figures.R        # publication figures (PDF + 600-dpi PNG/TIFF)
├── 05_tables.R         # result tables, descriptives, results_summary.md
└── run_all.R           # runs 01–05 and writes session info
Outputs/
├── Figures/            # fig1–fig7, figS1
├── Tables/             # model_results.{csv,xlsx}, model_diagnostics.csv,
│                       # descriptives.csv, data_checks.csv
├── results_summary.md  # ready-to-paste reporting sentences, caveats, figure legends
└── session_info.txt
```

`Outputs/Models/` holds the fitted model objects. It is regenerated on every run
and is not tracked by git.

## How to reproduce

Requires R ≥ 4.6 (developed with 4.6.1).

```r
# once, in the project root: installs the exact package versions from renv.lock
install.packages("renv")
renv::restore()
```

```sh
Rscript R/run_all.R
```

The full run takes roughly 15–25 minutes; the spatial refits in `03_spatial.R`
take most of that time. Each numbered script can also be run on its own,
because it reads its inputs from disk. All paths are relative (`here::here()`).
The random seed for the DHARMa simulations and figure jitter is set in
`00_setup.R`.

## Data

Each row is one mapping field (KFME cell, `POLE`). Columns are renamed to
English snake_case in `01_import_clean.R`:

| Original column | Clean name | Content |
|---|---|---|
| `ID`, `POLE` | `field_id`, `kfme_cell` | field ID and KFME grid code |
| `SMRK_m2` | `spruce_m2` → `spruce_km2` | remaining spruce forest area |
| `Pocet_nalezu` | `ndop_records` | number of NDOP (species occurrence database) records |
| `prob_1h_mean` | `prob_occ` | modelled occurrence probability 2014–2017 |
| `Plocha kurovcove kalamity celkem (m2)` | `calamity_m2` | bark-beetle calamity area (= dry + harvested) |
| `Plocha_souse (m2)` | `dry_m2` | standing dead spruce area |
| `Plocha_tezby (m2)` | `harvest_m2` | salvage-logged area |
| `Vyber_2025` | `selected_2025` | field selected for 2025 (yes/no) |
| `Mapovatel` | `observer` | surveyor(s) |
| `Datum_1`, `Datum_2` | `date_1`, `date_2` | visit dates (Excel serials in the source) |
| `Spontanne_v`, `Provokovano_v` | `spont_v`, `prov_v` | spontaneous / playback-provoked detections in visit *v* |
| `Termin_v` | `count_v` | visit total (= spontaneous + provoked) |
| `Celkem_zjisteno_oresniku` | `count_total` | **maximum** of the two visit totals |
| `Doba_2` | `time_of_day_2` | **time of day** of visit 2 (morning / afternoon-evening / both) – not provocation duration |
| `Zarizeni` | `device` | playback device (free text) |
| `Hlasitost` | `loudness_db` | device loudness at 1 m and maximum volume (dB; 63.5–119; 9 missing) |
| `DELKA_LINII` | `line_km` | total transect length walked (already in **km**; 1.7–15.4; 2 missing) |
| `X`, `Y` | `x_sjtsk`, `y_sjtsk` → `lon`, `lat` | field centroid |

### Coordinates

`X`/`Y` are **S-JTSK / Krovak East North (EPSG:5514)**. Both are negative;
`X` is the easting (more negative = further west) and `Y` the northing (more
negative = further south). `01_import_clean.R` converts them to WGS84 with
`sf`. It then checks every field against the centroid implied by its KFME
code: the maximum deviation is < 1 m. The pipeline stops if this check ever
fails, for example because of an axis swap in a future export. Latitude and
longitude in WGS84 are used for H10/H11, because the Krovak grid is rotated
against true north.

### Count unit

In every field `count_total = max(count_1, count_2)`, and
`count_v = spont_v + prov_v`. The source file does not say whether these are
**individuals or territories**. Figures and text currently say "territories";
to change this, edit `count_unit` in `00_setup.R`.

### Data cleaning decisions

All flagged records are listed in `Outputs/Tables/data_checks.csv`.

- "NO DATA" and "–" are treated as missing.
- Visit 2 was not done in 10 fields, but `Termin_2` was entered as 0. It is
  set to NA, and these fields are excluded from total-count models.
- Field 6048dc: the device note says visit 1 was done without playback, so
  its provoked visit-1 count is set to NA.
- Disturbance proportions = area / (remaining spruce + calamity area). The
  calamity area exceeds the remaining spruce in 8 fields, which the previous
  analysis dropped.
- Date models use only visits in the protocol windows (visit 1:
  February–April, visit 2: June–July). Two fields have visit 1 dated after
  visit 2.

## Effort correction

Transect length varied widely. 77 of 386 fields fall outside the 4–6 km
protocol (22 below, 55 above), 2 have no length, and none has zero length.
Every count model therefore uses the **raw count with
`offset(log(line_km))`**, so coefficients act on **density (count per km of
transect)** while the response keeps its count distribution. Effects are
reported as incidence rate ratios (IRR) and plotted per 1 km of transect.

Models are fitted with `glmmTMB`:

1. Poisson and negative binomial (NB2, the same likelihood as `MASS::glm.nb`)
   are fitted and compared by AIC and DHARMa dispersion tests.
2. NB2 is used unless its dispersion parameter runs to the Poisson boundary.
3. A zero-inflated NB is fitted only if DHARMa detects zero-inflation in the
   NB residuals, and kept only if it lowers AIC by more than 2.

Presence/absence models are binomial GLMs with a **cloglog** link and the same
offset. The detection rate then scales with transect length, and exp(β) is a
rate ratio.

## Hypotheses

| | Hypothesis | Response | Model |
|---|---|---|---|
| H1 | More spruce forest → more nutcrackers | total count / km | NB2 + offset |
| H2 | More NDOP records → more | total count / km | NB2 + offset |
| H3 | Higher modelled occurrence probability → more | total count / km | NB2 + offset |
| H4 | Larger bark-beetle calamity proportion → fewer | total count / km | NB2 + offset |
| H5 | Larger standing-dead proportion → fewer | total count / km | NB2 + offset |
| H6 | Larger salvage-logged proportion → fewer | total count / km | NB2 + offset |
| H7 | Date of visit 1 affects detections | visit-1 count / km | NB2 + offset |
| H8 | Date of visit 2 affects detections (adjusted for time of day) | visit-2 count / km | NB2 + offset |
| H9 | Louder playback device → more provoked detections, per visit | provoked count / km | NB2/Poisson + offset + (1 \| observer) |
| H10 | Further north → more nutcrackers | presence; total count / km | cloglog + offset; NB2 + offset |
| H11 | Further west → more nutcrackers | presence; total count / km | cloglog + offset; NB2 + offset |

H10 and H11 are fitted separately and jointly; the table reports the
latitude–longitude correlation and VIF. For H9, the observer is used as a
random effect. Device names are free text, and almost every observer used a
single device, so a device effect cannot be separated from an observer effect.
H9 sensitivity runs exclude the 119-dB megaphone and drop the random effect.

## Spatial autocorrelation

For every model, `03_spatial.R` tests the residuals with:

- DHARMa `testSpatialAutocorrelation()` (inverse-distance weights);
- Moran's I (`spdep`) on a symmetric 8-nearest-neighbour graph;
- a Moran's I correlogram in 10-km distance classes (Figure S1).

If any of them gives p < 0.05 (for the correlogram, the first 0–10 km class,
which catches dependence between adjacent fields that the global tests can
miss), the model is refitted with an **exponential
spatial covariance** on the field centroids, `exp(pos + 0 | group)` in
glmmTMB, and the tests are repeated. The fixed effects are unchanged.
Latitude and longitude stay as linear fixed effects; no spatial smooth of the
coordinates is used, because that would absorb the gradients H10/H11 test.
Only residual dependence is modelled. Non-spatial and spatial estimates are
both reported (`fit_type` in `model_results.csv`) so that spatial confounding
can be judged.

## Outputs

- `Outputs/Tables/model_results.csv` / `.xlsx`: one row per term, model and
  fit type. Columns cover estimate, SE, 95% CI, IRR or rate ratio, z, p, NB2
  θ, AIC, McFadden pseudo-R², spatial range, and spatial-autocorrelation tests
  before and after correction. `is_final = TRUE` marks the reported fit.
- `Outputs/Tables/model_diagnostics.csv`: DHARMa dispersion, zero-inflation,
  uniformity and outlier tests, plus Poisson/NB2/ZINB AIC.
- `Outputs/Tables/descriptives.csv`: sample sizes, effort distribution, raw
  and density summaries, X–Y correlation and VIF.
- `Outputs/results_summary.md`: one paragraph per hypothesis, caveats and
  figure legends.
- `Outputs/Figures/`: figures at 175 mm (double column) as PDF, 600-dpi PNG
  and TIFF.

## License

MIT – see [LICENSE](LICENSE).
