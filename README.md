# Nutcracker (*Nucifraga caryocatactes*) Data Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains scripts, data, and visualizations for analyzing [nutcracker](https://portal.nature.cz/w/druh-1281#/) survey data in relation to forest characteristics, disturbances, and observation conditions. The analyses focus on exploring hypotheses about factors affecting the number of observed ořešníks (*Celkem_zjisteno_oresniku*) using **negative binomial models (NB)** and visualizing predictions with **ggplot2**.

---

## Table of Contents

- [Contents](#contents)
- [Requirements](#requirements)
- [Data](#data)
- [Hypotheses](#hypotheses)
- [Analysis](#analysis)
- [Visualization](#visualization)
- [Outputs](#outputs)
- [How to Run](#how-to-run)
- [Notes](#notes)
- [License](#license)

---

## Contents
Data/  
├── Input/ # Raw Excel/CSV files  
└── Processed/ # Cleaned and processed CSV

Outputs/  
├── Analyses_nuc/ # Figures of predicted effects with observed points  
└── Tables/ # CSV/Excel tables of model coefficients

Scripts/  
└── analysis_nb.R # Main R script for negative binomial modeling

---
## Requirements

- R >= 4.3  
- Packages:

```r
tidyverse
MASS
ggplot2
patchwork
openxlsx
sjPlot
broom
report
```
Install missing packages with:
```r
install.packages(c("tidyverse","MASS","ggplot2","patchwork","openxlsx","sjPlot","broom","report"))
```

## Data
Oresnik_statistika.xlsx – main input data containing:
| Column                     | Description                            |
| -------------------------- | -------------------------------------- |
| `SMRK_m2`                  | Spruce forest area (m²)                |
| `Pocet_nalezu`             | Number of NDOP findings                |
| `prob_1h_mean`             | Modeled probability 2014–2017          |
| `Plocha_kurovcove_kal_m2`  | Area of bark beetle calamity           |
| `Plocha_souse_m2`          | Dead/dry wood area                     |
| `Plocha_tezby_m2`          | Harvested area                         |
| `Celkem_zjisteno_oresniku` | Total observed nutcrackers             |
| `Datum_1`, `Datum_2`       | Survey dates (DD/MM/YYYY)              |
| `Doba_2`                   | Provocation duration for second survey |

## Hypotheses
1) More spruce forest → more nutcrackers
2) More NDOP findings → more nutcrackers
3) Higher modeled probability → more nutcrackers
4) Higher proportion of calamity → fewer nutcrackers
5) Higher proportion of dead/dry trees → fewer nutcrackers
6) Higher proportion of harvested forest → fewer nutcrackers
7) Survey date in first term → effect on number of nutcrackers
8) Survey date in second term and provocation duration → effect on number of nutcrackers

## Analysis
- Negative binomial models are used to account for overdispersion.
- Missing data is ignored using na.exclude.
- Date predictors are converted to numeric days (as.numeric(Date)) for modeling.
- Interaction terms between forest area and disturbance variables can be included.
- Model coefficients are exported as CSV and readable tables.
- Visualization
- Predicted NB model effects are plotted with ggplot2.
- Plots include:
  - Observed data points (black, semi-transparent)
  - Predicted curve (blue line)
  - 95% confidence interval (shaded blue ribbon)
- Multiple hypotheses can be combined into panels using patchwork.
