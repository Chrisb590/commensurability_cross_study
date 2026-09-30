# Run the full analysis
#
# Sources every script in R/ one after the other. Inputs are read from data/
# and all figures and tables are written to output/ (created if missing).
#
# Usage (from the repository root):
#   Rscript run_analysis.R
# or, from an R session:
#   source("run_analysis.R")

library(cli)
library(here)

scripts <- c(
  "beta_diversity.R",
  "linear_models_vs_linear_mixed_models.R",
  "taxonomic_resolution_analysis.R",
  "patefield_null_model_analysis.R",
  "vaznull_model_analysis.R",
  "interaction_recording_and_edge_weight_units.R",
  "sampling_intensity_and_network_topology.R",
  "rarefaction.R",
  "map_of_locations.R"
)

# Each script is sourced in its own environment so that objects defined in one
# script cannot leak into the next.
run_script <- function(script) {
  cli_progress_step("Running {.file R/{script}}")
  source(here("R", script), local = new.env(), encoding = "UTF-8")
}

# When run non-interactively (e.g. Rscript), send the plots that scripts print()
# to a null device instead of creating Rplots.pdf in the working directory.
if (!interactive()) {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
}

cli_h1("Commensurability and cross-study inference")
cli_alert_info("Outputs will be written to {.path {here('output')}}")

for (script in scripts) run_script(script)

cli_alert_success("All {length(scripts)} scripts completed")
