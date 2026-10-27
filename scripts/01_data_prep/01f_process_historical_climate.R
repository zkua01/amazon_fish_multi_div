# =============================================================================
# Script:  01f_process_historical_climate.R
# Project: Amazon multidimensional biodiversity
# Purpose: Extract paleoclimate rasters (CCSM, MIROC, MPI-ESM-P) per
#          subdrainage basin, compute climate anomalies (current minus past),
#          average across models, and reduce via PCA.
#
# Inputs:
#   output/intermediate/basins_climate_current.rds   (from 01e)
#   data/raw/env/historical/past_CCSM/*.tif
#   data/raw/env/historical/past_MIROC/*.tif
#   data/raw/env/historical/past_MPI/*.tif
#
# Outputs:
#   output/intermediate/basins_historical_climate.rds
#   output/intermediate/pca_currpast.csv
#
# Note: Paleoclimate rasters from Varela et al. (2015) PLOS ONE.
#       See data/raw/README_download_instructions.md
#
# Dependencies: terra, sf, exactextractr, tidyverse, paran, here
# =============================================================================

library(terra)
library(sf)
library(exactextractr)
library(tidyverse)
library(paran)
library(here)

set.seed(123)

# --- Paths -------------------------------------------------------------------
inter_dir    <- here("output", "intermediate")
historical_dir <- here("data", "raw", "env", "historical")

# --- Helper function -----------------------------------------------------------
read_transform_rast <- function(path, target_crs = "EPSG:5880", method = "bilinear") {
  r <- rast(path)
  if (!is.na(crs(r)) && !grepl("5880", crs(r))) {
    message("Reprojecting raster: ", basename(path))
    r <- project(r, target_crs, method = method)
  }
  r
}

# --- 1. Load basins with current climate (from 01e) -----------------------------
message("Loading basins with current climate data...")
basins <- readRDS(file.path(inter_dir, "basins_climate_current.rds"))

# --- 2. Extract CCSM paleoclimate ------------------------------------------------
message("Extracting CCSM paleoclimate rasters...")
ccsm_dir <- file.path(historical_dir, "past_CCSM")

ccsm_vars <- c(
  cctemp = "cctemp.tif",
  cctmax = "cctmax.tif",
  cctmin = "cctmin.tif",
  ccprep = "ccprep.tif",
  ccpmax = "ccpmax.tif",
  ccpmin = "ccpmin.tif"
)

for (var in names(ccsm_vars)) {
  r <- read_transform_rast(file.path(ccsm_dir, ccsm_vars[var]))
  basins[[var]] <- exact_extract(r, basins, fun = "mean")
  rm(r)
}
gc()

# --- 3. Extract MIROC paleoclimate ------------------------------------------------
message("Extracting MIROC paleoclimate rasters...")
miroc_dir <- file.path(historical_dir, "past_MIROC")

miroc_vars <- c(
  mrtemp = "mrtemp.tif",
  mrtmax = "mrtmax.tif",
  mrtmin = "mrtmin.tif",
  mrprep = "mrprep.tif",
  mrpmax = "mrpmax.tif",
  mrpmin = "mrpmin.tif"
)

for (var in names(miroc_vars)) {
  r <- read_transform_rast(file.path(miroc_dir, miroc_vars[var]))
  basins[[var]] <- exact_extract(r, basins, fun = "mean")
  rm(r)
}
gc()

# --- 4. Extract MPI-ESM-P paleoclimate ---------------------------------------------
message("Extracting MPI-ESM-P paleoclimate rasters...")
mpi_dir <- file.path(historical_dir, "past_MPI")

mpi_vars <- c(
  metemp = "metemp.tif",
  metmax = "metmax.tif",
  metmin = "metmin.tif",
  meprep = "meprep.tif",
  mepmax = "mepmax.tif",
  mepmin = "mepmin.tif"
)

for (var in names(mpi_vars)) {
  r <- read_transform_rast(file.path(mpi_dir, mpi_vars[var]))
  basins[[var]] <- exact_extract(r, basins, fun = "mean")
  rm(r)
}
gc()

# --- 5. Save intermediate ----------------------------------------------------------
saveRDS(basins, file.path(inter_dir, "basins_historical_climate.rds"))
message("Basins with historical climate saved.")

# --- 6. Compute climate anomalies (current - past) ----------------------------------
message("Computing climate anomalies...")
basins_his <- st_drop_geometry(basins)

basins_his <- basins_his |>
  mutate(
    # CCSM anomalies
    diff_cctemp = temp - cctemp,
    diff_cctmax = tmax - cctmax,
    diff_cctmin = tmin - cctmin,
    diff_ccprep = prep - ccprep,
    diff_ccpmax = pmax - ccpmax,
    diff_ccpmin = pmin - ccpmin,
    # MIROC anomalies
    diff_mrtemp = temp - mrtemp,
    diff_mrtmax = tmax - mrtmax,
    diff_mrtmin = tmin - mrtmin,
    diff_mrprep = prep - mrprep,
    diff_mrpmax = pmax - mrpmax,
    diff_mrpmin = pmin - mrpmin,
    # MPI anomalies
    diff_metemp = temp - metemp,
    diff_metmax = tmax - metmax,
    diff_metmin = tmin - metmin,
    diff_meprep = prep - meprep,
    diff_mepmax = pmax - mepmax,
    diff_mepmin = pmin - mepmin
  )

# --- 7. Average anomalies across the three models ------------------------------------
basins_his <- basins_his |>
  mutate(
    temp_change = (diff_mrtemp + diff_cctemp + diff_metemp) / 3,
    tmax_change = (diff_mrtmax + diff_cctmax + diff_metmax) / 3,
    tmin_change = (diff_mrtmin + diff_cctmin + diff_metmin) / 3,
    prep_change = (diff_mrprep + diff_ccprep + diff_meprep) / 3,
    pmax_change = (diff_mrpmax + diff_ccpmax + diff_mepmax) / 3,
    pmin_change = (diff_mrpmin + diff_ccpmin + diff_mepmin) / 3
  )

# --- 8. PCA on averaged climate anomalies ----------------------------------------------
message("Running PCA: Climate anomaly (CurrPast)...")
currpast <- basins_his |>
  dplyr::select(temp_change, tmax_change, tmin_change,
                prep_change, pmax_change, pmin_change)
currpast_sc <- as.data.frame(scale(currpast))

paran_currpast <- paran(as.matrix(currpast_sc), iterations = 5000)
pca_currpast   <- prcomp(currpast_sc)

pca_currpast_scores <- as.data.frame(pca_currpast$x) |>
  dplyr::select(1:paran_currpast$Retained) |>
  dplyr::rename_with(~ paste(.x, "CurrPast", sep = "_")) |>
  rownames_to_column(var = "SBD_CD")

write.csv(pca_currpast_scores, file.path(inter_dir, "pca_currpast.csv"), row.names = FALSE)
message("PCA CurrPast: ", paran_currpast$Retained, " components retained.")

# --- Done --------------------------------------------------------------------------------
message("01f complete. Outputs saved to: ", inter_dir)
