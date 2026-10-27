# =============================================================================
# Script:  01e_extract_climate_rasters.R
# Project: Amazon multidimensional biodiversity
# Purpose: Extract current-climate raster data (WorldClim, MODIS NPP, runoff)
#          per subdrainage basin, then reduce multicollinearity via PCA for
#          three variable groups: Water, Temperature, and Energy.
#
# Inputs:
#   data/raw/SubDrainageShapefile/SubDrainageShapefile.shp   (basin polygons)
#   data/raw/env/climate/wc2.1_10m_bio/*.tif                 (WorldClim Bio)
#   data/raw/env/climate/wc2.1_10m_srad/*.tif                (solar radiation)
#   data/raw/env/climate/wc2.1_10m_vapr/*.tif                (water vapour)
#   data/raw/env/climate/PET-Global-ET0_v3_monthly/*.tif     (PET)
#   data/raw/env/climate/aet_monthly/                        (AET, nested)
#   data/raw/env/climate/MOD17A3HGF_NPP_00_25.tif            (NPP)
#   data/raw/env/habitat/Runoff_Sam/cmp_sam_ro.grd           (runoff)
#
# Outputs:
#   output/intermediate/basins_climate_current.rds  (full extracted table)
#   output/intermediate/pca_water.csv
#   output/intermediate/pca_temp.csv
#   output/intermediate/pca_energy.csv
#
# Note: Large rasters must be downloaded separately.
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
basin_path  <- here("data", "raw", "SubDrainageShapefile", "SubDrainageShapefile.shp")
climate_dir <- here("data", "raw", "env", "climate")
habitat_dir <- here("data", "raw", "env", "habitat", "Runoff_Sam")
inter_dir   <- here("output", "intermediate")
dir.create(inter_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("temp_swap", showWarnings = FALSE)

terraOptions(memfrac = 0.8, tempdir = "temp_swap", verbose = TRUE)

# --- Helper functions ----------------------------------------------------------

# Read and reproject a vector file
read_transform <- function(path, target_crs = 5880) {
  obj <- st_read(path, quiet = TRUE)
  if (st_crs(obj) != target_crs) {
    obj <- st_transform(obj, target_crs)
  }
  obj
}

# Read and reproject a raster file
read_transform_rast <- function(path, target_crs = "EPSG:5880", method = "bilinear") {
  r <- rast(path)
  if (!is.na(crs(r)) && !grepl("5880", crs(r))) {
    message("Reprojecting raster: ", basename(path))
    r <- project(r, target_crs, method = method)
  }
  r
}

# Calculate summary statistics from a monthly data frame
calc_stats <- function(monthly_df) {
  data.frame(
    mean   = rowMeans(monthly_df, na.rm = TRUE),
    annual = rowSums(monthly_df, na.rm = TRUE),
    max    = apply(monthly_df, 1, max, na.rm = TRUE),
    min    = apply(monthly_df, 1, min, na.rm = TRUE),
    cv     = apply(monthly_df, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  )
}

# Process monthly climate rasters and extract basin-level statistics
process_climate_data <- function(variable_name, file_path, basins_object,
                                 months, file_structure = "flat") {
  monthly_results <- data.frame(matrix(nrow = nrow(basins_object), ncol = 0))
  basins_extent   <- terra::ext(basins_raw) + 0.1
  
  for (month in months) {
    message("Processing ", variable_name, " month ", month, "...")
    
    if (file_structure == "nested") {
      folder_path <- file.path(file_path, paste0(variable_name, "_", as.numeric(month)))
      file        <- file.path(folder_path, "w001001.adf")
    } else {
      file <- file.path(file_path, paste0(variable_name, "_", month, ".tif"))
    }
    message("Reading: ", file)
    
    if (!file.exists(file)) stop("File not found: ", file)
    
    r <- crop(rast(file), basins_extent)
    if (!is.na(crs(r)) && !grepl("EPSG:5880", crs(r))) {
      message("Reprojecting cropped raster: ", basename(file))
      r <- project(r, "EPSG:5880", method = "bilinear")
    }
    
    message("Extracting values...")
    month_values <- exact_extract(r, basins_object, fun = "mean")
    monthly_results[[paste0("month_", month)]] <- month_values
    
    rm(r, month_values)
    gc()
    message("Completed month ", month)
  }
  
  calc_stats(monthly_results)
}


# --- 1. Load basin shapefile ---------------------------------------------------
message("Loading basin shapefile...")
basins_raw <- st_read(basin_path)
basins <- read_transform(basin_path)
basins <- basins[order(basins$SBD_CD), ]
rownames(basins) <- basins$SBD_CD

# --- 2. WorldClim Bio variables (current climate) ------------------------------
message("Extracting WorldClim Bio variables...")

wc_files <- c(
  temp = "wc2.1_10m_bio_1.tif",
  tmax = "wc2.1_10m_bio_5.tif",
  tmin = "wc2.1_10m_bio_6.tif",
  tcv  = "wc2.1_10m_bio_4.tif",
  prep = "wc2.1_10m_bio_12.tif",
  pmax = "wc2.1_10m_bio_13.tif",
  pmin = "wc2.1_10m_bio_14.tif",
  pcv  = "wc2.1_10m_bio_15.tif"
)

bio_dir <- file.path(climate_dir, "wc2.1_10m_bio")
for (var in names(wc_files)) {
  r <- read_transform_rast(file.path(bio_dir, wc_files[var]))
  basins[[var]] <- exact_extract(r, basins, fun = "mean")
  rm(r)
}
gc()

# --- 3. Monthly climate variables (srad, wvp, PET, AET) ------------------------
message("Processing monthly climate variables...")
months <- sprintf("%02d", 1:12)

stats_srad <- process_climate_data("wc2.1_10m_srad",
                                   file.path(climate_dir, "wc2.1_10m_srad"),
                                   basins, months)

stats_wvp <- process_climate_data("wc2.1_10m_vapr",
                                  file.path(climate_dir, "wc2.1_10m_vapr"),
                                  basins, months)

stats_pet <- process_climate_data("et0_v3",
                                  file.path(climate_dir, "PET-Global-ET0_v3_monthly"),
                                  basins, months)

stats_aet <- process_climate_data("aet",
                                  file.path(climate_dir, "aet_monthly"),
                                  basins, months,
                                  file_structure = "nested")

# Add monthly-derived statistics to basins
basins <- basins |>
  mutate(
    srad_mean  = stats_srad$mean,
    srad_cv    = stats_srad$cv,
    wvp_mean   = stats_wvp$mean,
    wvp_cv     = stats_wvp$cv,
    aet_annual = stats_aet$annual,
    aet_max    = stats_aet$max,
    aet_min    = stats_aet$min,
    aet_cv     = stats_aet$cv,
    pet_annual = stats_pet$annual,
    pet_max    = stats_pet$max,
    pet_min    = stats_pet$min,
    pet_cv     = stats_pet$cv
  )

# --- 4. MODIS NPP ---------------------------------------------------------------
message("Extracting MODIS NPP...")
mod17_npp <- read_transform_rast(file.path(climate_dir, "MOD17A3HGF_NPP_00_25.tif"))
basins$npp_mean <- exact_extract(mod17_npp, basins, fun = "mean") |> rowMeans()
basins$npp_sd   <- exact_extract(mod17_npp, basins,
                                 fun = function(v, cf) apply(v, 2, sd, na.rm = TRUE)) |>
  rowMeans(na.rm = TRUE)
basins$npp_cv <- basins$npp_sd / basins$npp_mean
rm(mod17_npp)
gc()

# --- 5. Runoff -------------------------------------------------------------------
message("Extracting runoff...")
runoff_raster <- read_transform_rast(file.path(habitat_dir, "cmp_sam_ro.grd"))
basins$runoff <- exact_extract(runoff_raster, basins, fun = "mean")
rm(runoff_raster)
gc()

# --- 6. Save intermediate extracted table -----------------------------------------
saveRDS(basins, file.path(inter_dir, "basins_climate_current.rds"))
message("Basins climate table saved.")

# --- 7. PCA: Water (precipitation + runoff) --------------------------------------
message("Running PCA: Water...")
basins_curr <- st_drop_geometry(basins)

water_vars <- basins_curr |> dplyr::select(prep, pmax, pmin, pcv, runoff)
water_sc   <- as.data.frame(scale(water_vars))

paran_water <- paran(as.matrix(water_sc), iterations = 5000)
pca_water   <- prcomp(water_sc)

pca_water_scores <- as.data.frame(pca_water$x) |>
  dplyr::select(1:paran_water$Retained) |>
  dplyr::rename_with(~ paste(.x, "Water", sep = "_")) |>
  rownames_to_column(var = "SBD_CD")

write.csv(pca_water_scores, file.path(inter_dir, "pca_water.csv"), row.names = FALSE)
message("PCA Water: ", paran_water$Retained, " components retained.")

# --- 8. PCA: Temperature (temp + solar radiation) ---------------------------------
message("Running PCA: Temperature...")
temp_vars <- basins_curr |> dplyr::select(temp, tmax, tmin, tcv, srad_mean, srad_cv)
temp_sc   <- as.data.frame(scale(temp_vars))

paran_temp <- paran(as.matrix(temp_sc), iterations = 5000)
pca_temp   <- prcomp(temp_sc)

pca_temp_scores <- as.data.frame(pca_temp$x) |>
  dplyr::select(1:paran_temp$Retained) |>
  dplyr::rename_with(~ paste(.x, "Temp", sep = "_")) |>
  rownames_to_column(var = "SBD_CD")

write.csv(pca_temp_scores, file.path(inter_dir, "pca_temp.csv"), row.names = FALSE)
message("PCA Temperature: ", paran_temp$Retained, " components retained.")

# --- 9. PCA: Energy (AET + PET + NPP + WVP) --------------------------------------
message("Running PCA: Energy...")
energy_vars <- basins_curr |>
  dplyr::select(aet_annual, aet_max, aet_min, aet_cv,
                pet_annual, pet_max, pet_min, pet_cv,
                npp_mean, npp_sd, npp_cv)
energy_sc <- as.data.frame(scale(energy_vars))

paran_energy <- paran(as.matrix(energy_sc), iterations = 5000)
pca_energy   <- prcomp(energy_sc)

pca_energy_scores <- as.data.frame(pca_energy$x) |>
  dplyr::select(1:paran_energy$Retained) |>
  dplyr::rename_with(~ paste(.x, "Energy", sep = "_")) |>
  rownames_to_column(var = "SBD_CD")

write.csv(pca_energy_scores, file.path(inter_dir, "pca_energy.csv"), row.names = FALSE)
message("PCA Energy: ", paran_energy$Retained, " components retained.")

# --- Done --------------------------------------------------------------------------
message("01e complete. Outputs saved to: ", inter_dir)
