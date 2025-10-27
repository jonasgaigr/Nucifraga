#----------------------------------------------------------#
# Load packages -----
#----------------------------------------------------------#
packages <- c(
  "tidyverse",
  "janitor",
  "sf", 
  "sp", 
  "proj4", 
  "openxlsx",
  "fuzzyjoin", 
  "remotes",
  "ggtext"
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
raw <- openxlsx::read.xlsx("Data/Input/Oresnik_statistika.xlsx", sheet = 1, colNames = TRUE)

df <- raw %>%
  rename_with(~ str_replace_all(., "\\.", "_")) %>%     # replace "." with "_"
  rename_with(~ str_remove_all(., "[()]"))     


# Write processed data ----
write_csv2(
  df,
  "Data/Processed/data_nuc.csv"
)

