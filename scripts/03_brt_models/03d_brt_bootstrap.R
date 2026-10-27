# =============================================================================
# Script:  03d_brt_bootstrap.R
# Project: Amazon multidimensional biodiversity
# Purpose: Fit final BRT models with optimized hyperparameters and quantify
#          uncertainty via spatial block bootstrap (1000 iterations for
#          relative influence and R²; 200 for PDPs and Friedman's H).
#
# Inputs:
#   data/processed/dataprepared_filtered.rds   (from 03b)
#   output/intermediate/brt_best_params.csv    (from 03c)
#   data/raw/SubDrainageShapefile/SubDrainageShapefile.shp
#
# Outputs:
#   output/intermediate/brt_final_results.rds
#   output/intermediate/brt_influence_td.csv
#   output/intermediate/brt_influence_pd.csv
#   output/intermediate/brt_influence_fd.csv
#   output/intermediate/brt_influence_md.csv
#
# Dependencies: gbm3, dismo, blockCV, sf, tidyverse, Rmisc, here
# =============================================================================

library(gbm3)
library(dismo)
library(blockCV)
library(sf)
library(tidyverse)
library(Rmisc)
library(here)

source(here("scripts", "utils", "brt_function.R"))
source(here("scripts", "utils", "brt_helper_fn.R"))

# --- Paths ---
input_dir  <- here("data", "processed")
inter_dir  <- here("output", "intermediate")
basin_path <- here("data", "raw", "SubDrainageShapefile", "SubDrainageShapefile.shp")

# --- 1. Load data ---
dataprepared_filt <- readRDS(file.path(input_dir, "dataprepared_filtered.rds"))
best_params <- read.csv(file.path(inter_dir, "brt_best_params.csv"))
basins <- st_read(basin_path, quiet = TRUE) |>
  st_transform(crs = 4326) |> arrange(SBD_CD)

responses <- c("muDiv", "logSR", "sesFD", "sesPD")
dp <- dataprepared_filt |> na.omit() |>
  left_join(dplyr::select(basins, SBD_CD, geometry), by = "SBD_CD")

div <- dplyr::select(dp, all_of(responses))
polygon <- dplyr::select(dp, SBD_CD, geometry) |> st_as_sf() |> st_make_valid()
coord <- dplyr::select(dp, CENT_X, CENT_Y)
env <- dplyr::select(dp, AREA:sample_area) |>
  mutate(across(AREA, log1p)) |>
  mutate(across(everything(), scale)) |>
  mutate(across(everything(), as.numeric))

df <- bind_cols(div, env, coord)
predictors <- names(env)

# --- 2. Run analysis per response ---
create_spatial_folds <- function(data, cv, polygon) {
  data_sf <- st_as_sf(data, coords = c("CENT_X", "CENT_Y"), crs = 4326)
  cv_spatial(x = data_sf, k = cv, user_blocks = polygon,
             selection = "random", iteration = 100, progress = FALSE)$folds_list
}

run_analysis <- function(data, response, predictors, params_df, polygon,
                         n_influence_reps = 1000, n_pdp_reps = 200, n_h_reps = 200) {
  params <- params_df |> filter(div == response)
  message("\n=== Processing: ", response, " ===")
  message("  n_trees=", params$n_trees, " lr=", params$lr,
          " tc=", params$tc, " bf=", params$bf)
  
  spatial_folds <- create_spatial_folds(data, cv = 10, polygon = polygon)
  
  # --- Relative influence (1000 bootstrap) ---
  message("  Bootstrapping relative influence...")
  set.seed(324)
  influence_results <- manual_spatial_bootstrap(
    data, response, predictors, params, spatial_folds,
    calc_influence, R = n_influence_reps
  )
  influence_ci <- t(apply(influence_results, 2, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) < 10) return(c(NA, NA, NA))
    Rmisc::CI(x_clean)
  }))
  colnames(influence_ci) <- c("lower", "mean", "upper")
  rownames(influence_ci) <- predictors
  
  # --- R² (1000 bootstrap) ---
  message("  Bootstrapping R²...")
  set.seed(324)
  r2_results <- manual_spatial_bootstrap(
    data, response, predictors, params, spatial_folds,
    calc_r2, R = n_influence_reps
  )
  r2_clean <- r2_results[, 1][!is.na(r2_results[, 1])]
  r2_ci <- if (length(r2_clean) < 10) c(NA, NA, NA) else Rmisc::CI(r2_clean)
  names(r2_ci) <- c("lower", "mean", "upper")
  
  # --- PDPs (200 bootstrap) ---
  message("  Calculating PDPs...")
  pdp_list <- list()
  pdp_results_list <- list()
  for (var in predictors) {
    set.seed(324)
    pdp_results <- manual_spatial_bootstrap(
      data, response, predictors, params, spatial_folds,
      calc_pdp, R = n_pdp_reps,
      additional_args = list(focal_var = var, grid_points = 50)
    )
    pdp_results_list[[var]] <- pdp_results
    
    var_range <- range(data[[var]], na.rm = TRUE)
    grid_vals <- if (diff(var_range) == 0) rep(var_range[1], 50)
    else seq(var_range[1], var_range[2], length.out = 50)
    
    pdp_list[[var]] <- data.frame(
      variable = var, value = grid_vals,
      mean  = apply(pdp_results, 2, mean, na.rm = TRUE),
      lower = apply(pdp_results, 2, quantile, 0.025, na.rm = TRUE),
      upper = apply(pdp_results, 2, quantile, 0.975, na.rm = TRUE)
    )
  }
  
  # --- Friedman's H (200 bootstrap) ---
  message("  Bootstrapping Friedman's H...")
  set.seed(324)
  h_results <- manual_spatial_bootstrap(
    data, response, predictors, params, spatial_folds,
    calc_friedman_h, R = n_h_reps
  )
  pairs <- combn(predictors, 2, simplify = FALSE)
  pair_names <- vapply(pairs, function(p) paste(p, collapse = " x "), character(1))
  h_ci <- t(apply(h_results, 2, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) < 10) return(c(NA, NA, NA))
    Rmisc::CI(x_clean)
  }))
  colnames(h_ci) <- c("lower", "mean", "upper")
  rownames(h_ci) <- pair_names
  
  list(
    response = response, parameters = params,
    influence = influence_ci, r2 = r2_ci,
    pdp = pdp_list, friedman_h = h_ci,
    bootstrap_raw = list(influence = influence_results, r2 = r2_results,
                         pdp = pdp_results_list, friedman_h = h_results)
  )
}

# --- 3. Execute ---
final_results <- list()
for (resp in responses) {
  tryCatch({
    final_results[[resp]] <- run_analysis(
      df, resp, predictors, best_params, polygon,
      n_influence_reps = 1000, n_pdp_reps = 200, n_h_reps = 200
    )
  }, error = function(e) {
    message("ERROR processing ", resp, ": ", e$message)
  })
}

# --- 4. Save ---
saveRDS(final_results, file.path(inter_dir, "brt_final_results.rds"))

for (resp in names(final_results)) {
  suffix <- case_when(
    resp == "logSR" ~ "td", resp == "sesPD" ~ "pd",
    resp == "sesFD" ~ "fd", resp == "muDiv" ~ "md"
  )
  final_results[[resp]]$influence |>
    as.data.frame() |>
    rownames_to_column(var = "predictor") |>
    arrange(desc(mean)) |>
    write.csv(file.path(inter_dir, paste0("brt_influence_", suffix, ".csv")),
              row.names = FALSE)
}

message("03d complete.")
