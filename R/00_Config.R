#----------------------------------------------------------#
# Load packages -----
#----------------------------------------------------------#
packages <- c(
  "tidyverse",
  "broom",
  "janitor",
  "sf", 
  "sp", 
  "proj4", 
  "openxlsx",
  "fuzzyjoin", 
  "remotes",
  "ggtext",
  "vegan",
  "ggplot2",
  "ggrepel",
  "ggforce",
  "MASS",
  "DHARMa",
  "patchwork",
  "sjPlot",
  "report"
)

# Standard package
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

#----------------------------------------------------------#
# Load data -----
#----------------------------------------------------------#
raw <- openxlsx::read.xlsx(
  "Data/Input/Oresnik_statistika.xlsx", 
  sheet = 1, 
  colNames = TRUE,
  detectDates = TRUE
  )

df <- raw %>%
  rename_with(~ str_replace_all(., "\\.", "_")) %>%     # replace "." with "_"
  rename_with(~ str_remove_all(., "[()]")) %>%
  # Replace Inf / -Inf with NA
  mutate(across(everything(), ~ ifelse(is.infinite(.x), NA, .x)))


# Write processed data ----
write_csv2(
  df,
  "Data/Processed/data_nuc.csv"
)

