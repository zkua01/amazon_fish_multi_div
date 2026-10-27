# =============================================================================
# Script:  03a_assemble_modelling_data.R
# Project: Amazon multidimensional biodiversity
# Purpose: Assemble the final environmental predictor table from PCA scores
#          and habitat variables, then join with the diversity metrics table
#          to produce the modelling-ready dataset consumed by BRT.
#
# Inputs:
#   output/intermediate/pca_water.csv             (from 01e)
#   output/intermediate/pca_temp.csv              (from 01e)
#   output/intermediate/pca_energy.csv            (from 01e)
#   output/intermediate/pca_currpast.csv          (from 01f)
#   output/intermediate/pca_elev.csv              (from 01g)
#   output/intermediate/habitat_position_vars.csv (from 01g)
#   data/processed/diversity_metrics_per_basin.csv (from 02e)
#
# Outputs:
#   data/processed/environmental_predictors.csv   (env table only, for reference)
#   data/processed/dataprepared.rds               (joined responses + predictors)
#
# Dependencies: tidyverse, here
# =============================================================================

library(tidyverse)
library(here)

# --- Paths ---
inter_dir   <- here("output", "intermediate")
input_dir   <- here("data", "processed")
output_dir  <- here("data", "processed")

# --- 1. Load PCA scores from earlier steps ---
pca_water    <- read.csv(file.path(inter_dir, "pca_water.csv"))
pca_temp     <- read.csv(file.path(inter_dir, "pca_temp.csv"))
pca_energy   <- read.csv(file.path(inter_dir, "pca_energy.csv"))
pca_currpast <- read.csv(file.path(inter_dir, "pca_currpast.csv"))
pca_elev     <- read.csv(file.path(inter_dir, "pca_elev.csv"))

# --- 2. Load habitat/position variables (from 01g) ---
habitat_vars <- read.csv(file.path(inter_dir, "habitat_position_vars.csv"))

# --- 3. Assemble environmental predictor table ---
envdata <- dplyr::select(habitat_vars, SBD_CD, CENT_X, CENT_Y, AREA) |>
  left_join(pca_water, by = "SBD_CD") |>
  left_join(pca_temp, by = "SBD_CD") |>
  left_join(pca_energy, by = "SBD_CD") |>
  left_join(pca_currpast, by = "SBD_CD") |>
  left_join(pca_elev, by = "SBD_CD") |>
  left_join(dplyr::select(habitat_vars, SBD_CD, PebasLake:sample_area), by = "SBD_CD")

write.csv(envdata,
          file.path(output_dir, "environmental_predictors.csv"),
          row.names = FALSE)
message("Environmental predictor table saved: ",
        nrow(envdata), " basins x ", ncol(envdata), " variables")

# --- 4. Load diversity metrics (from 02e) ---
div_tbl <- read.csv(file.path(input_dir, "diversity_metrics_per_basin.csv"))

# --- 5. Join diversity + environment → modelling dataset ---
dataprepared <- left_join(div_tbl, envdata, by = "SBD_CD")

# Validate: no unmatched rows
n_div <- nrow(div_tbl)
n_joined <- nrow(dataprepared)
stopifnot(n_div == n_joined)
message("Joined dataset: ", n_joined, " basins x ", ncol(dataprepared), " columns")

# --- 6. Save ---
saveRDS(dataprepared, file.path(output_dir, "dataprepared.rds"))
message("Modelling dataset saved: data/processed/dataprepared.rds")
