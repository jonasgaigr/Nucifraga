#----------------------------------------------------------#
# Run the full pipeline from the raw Excel file -----
#----------------------------------------------------------#
# Usage (from the project root, after `renv::restore()`):
#   Rscript R/run_all.R
# Runtime is roughly 25-40 minutes; the spatial refits in 03 dominate.

scripts <- c(
  "01_import_clean.R",
  "02_models.R",
  "03_spatial.R",
  "04_figures.R",
  "05_tables.R",
  "06_exploratory_pca.R"
)

for (script in scripts) {
  message("\n==== ", script, " ====")
  # A fresh environment per script, so each one only relies on what it reads
  # from disk and on 00_setup.R
  source(here::here("R", script), local = new.env())
}

source(here::here("R", "00_setup.R"))
writeLines(capture.output(sessionInfo()), here("Outputs", "session_info.txt"))
message("\nPipeline finished; session info in Outputs/session_info.txt")
