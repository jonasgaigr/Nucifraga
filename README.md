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
├── 06_exploratory_pca.R # exploratory PCA of the habitat variables
└── run_all.R           # runs 01–06 and writes session info
Outputs/
├── Figures/            # fig1–fig7, figS1–figS3
├── Tables/             # model_results.{csv,xlsx}, model_diagnostics.csv,
│                       # descriptives.csv, data_checks.csv,
│                       # selection_sensitivity.csv, pca_*.csv
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

The full run takes roughly 25–40 minutes; the spatial refits in `03_spatial.R`
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
| `SMRK_m2` | `spruce_m2` → `spruce_km2` | spruce forest area remaining after the calamity |
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
| `DELKA_LINII` | `line_km` | transect length walked per visit (already in **km**; 1.7–15.4; 2 missing) |
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
`count_v = spont_v + prov_v`. The counts are **territories** (confirmed by the
data owner; the source file itself does not say). Labels follow `count_unit`
in `00_setup.R`.

### Data cleaning decisions

All flagged records are listed in `Outputs/Tables/data_checks.csv`.

- "NO DATA" and "–" are treated as missing.
- Visit 2 was not done in 10 fields, but `Termin_2` was entered as 0. It is
  set to NA, and these fields are excluded from total-count models.
- Field 6048dc: the device note says visit 1 was done without playback, so
  its provoked visit-1 count is set to NA.
- Fields 5253ab and 5254bb: visit 1 was entered in August, after visit 2
  (confirmed entry errors). The month is corrected from 08 to 03 in
  `01_import_clean.R`. The corrected dates (30 and 23 March) fall within the
  same observer's other visit-1 dates (22 March–3 April), on days on which
  neighbouring fields were surveyed. The script stops if the source values
  change, and if any visit 1 is still dated after visit 2.
- Transect lengths > 10 km (13 fields) were confirmed as single-visit lengths,
  not the sum of both visits, and are kept.
- Device identities that one observer marked as uncertain (16 fields) were
  verified by the data owner and are kept.
- Disturbance proportions = area / (remaining spruce + calamity area);
  `SMRK_m2` is the spruce remaining after the calamity (confirmed). The
  calamity area exceeds the remaining spruce in 8 fields, which the previous
  analysis dropped.
- Date models use only visits in the protocol windows (visit 1:
  February–April, visit 2: June–July).

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

## Sampling design (`Vyber_2025`)

Fields selected for 2025 have about twice the mean count of the other fields.
Selection was most likely based on expected occurrence, so the habitat effects
may be partly a product of the design. The main models leave selection status
out; `02_models.R` adds it as a covariate in sensitivity models (`SEL` for the
selection effect alone, `H1_sel`–`H6_sel` and `H10_H11_dens_sel`). The
adjusted and unadjusted effects are compared in
`Outputs/Tables/selection_sensitivity.csv`, `results_summary.md` and
Figure S2. Because H2 (NDOP records) and H3 (modelled occurrence probability)
are themselves measures of expected occurrence, a weaker effect after
adjustment there may mean shared variation with the selection criteria, not
only a design artefact.

## Exploratory PCA

`06_exploratory_pca.R` replaces the earlier vegan PCA, which mixed counts with
predictors, used raw m² areas and grouped points by field ID. The new version
is a correlation-matrix PCA of the habitat predictors only, on the fields in
the total-count models:

- spruce area;
- NDOP records, log(1 + x);
- occurrence probability, log;
- standing dead and salvage-logged proportions, logit.

Calamity is left out because it equals dry + harvested area. Density (square
root), latitude and longitude are projected passively, so they do not shape
the axes. Components are retained by the broken-stick rule. An NB2/Poisson
model of the total count on PC1 and PC2, with the same offset and a spatial
refit, shows how much density follows the joint habitat gradient. Outputs:
`pca_variance.csv`, `pca_loadings.csv`, `pca_density_model.csv` and Figure S3.
This analysis is descriptive; the hypothesis tests are in 02–03.

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
  and density summaries, X–Y correlation and VIF, and summaries by selection
  status.
- `Outputs/Tables/selection_sensitivity.csv`: habitat and geographic effects
  with and without selection status as a covariate.
- `Outputs/Tables/pca_*.csv`: exploratory PCA (variance, loadings, density
  model).
- `Outputs/results_summary.md`: one paragraph per hypothesis, caveats and
  figure legends.
- `Outputs/Figures/`: figures as PDF, 600-dpi PNG and TIFF, mostly at 175 mm
  (double column).

## License

MIT – see [LICENSE](LICENSE).
