# =============================================================================
# Script:  03c_brt_optimize.R
# Project: Amazon multidimensional biodiversity
# Purpose: Optimize BRT hyperparameters (learning rate, tree complexity, bag
#          fraction) via grid search with 10-fold spatial cross-validation.
#          Spatial folds are defined by subdrainage basin polygons (blockCV)
#          to account for spatial autocorrelation.
#
# Inputs:
#   data/processed/dataprepared_filtered.rds  (from 03b)
#   data/raw/SubDrainageShapefile/SubDrainageShapefile.shp
#
# Outputs:
#   output/intermediate/brt_best_params.csv
#   output/intermediate/brt_optimization_results.rds
#
# Dependencies: gbm3, blockCV, sf, tidyverse, here
# =============================================================================

library(gbm3)
library(blockCV)
library(sf)
library(tidyverse)
library(here)

# --- Paths ---
input_dir   <- here("data", "processed")
basin_path  <- here("data", "raw", "SubDrainageShapefile", "SubDrainageShapefile.shp")
output_dir  <- here("output", "intermediate")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load data ---
message("Loading prepared data...")
dataprepared_filt <- readRDS(file.path(input_dir, "dataprepared_filtered.rds"))
basins <- st_read(basin_path, quiet = TRUE) |>
  st_transform(crs = 4326) |>
  arrange(SBD_CD)

# --- 2. Assemble modelling dataframe ---
dp <- dataprepared_filt |>
  na.omit() |>
  left_join(dplyr::select(basins, SBD_CD, geometry), by = "SBD_CD")

responses <- c("muDiv", "logSR", "sesFD", "sesPD")

div <- dplyr::select(dp, all_of(responses))
polygon <- dplyr::select(dp, SBD_CD, geometry) |> st_as_sf() |> st_make_valid()
coord <- dplyr::select(dp, CENT_X, CENT_Y)

# Predictors: log-transform area, then standardise all
env <- dplyr::select(dp, AREA:sample_area) |>
  mutate(across(AREA, log1p)) |>
  mutate(across(everything(), scale)) |>
  mutate(across(everything(), as.numeric))

df <- bind_cols(div, env, coord)
predictors <- names(env)

message("Modelling data: ", nrow(df), " basins, ",
        length(predictors), " predictors")

# --- 3. Define hyperparameter grid ---
hyper_grid <- expand.grid(
  lr = c(0.01, 0.005, 0.001),
  bf = seq(0.5, 0.8, 0.15),
  tc = seq(2, 5, 1)
)

# --- 4. Spatial fold creation (once per response) ---
create_spatial_folds <- function(data, cv, polygon,
                                 coords_columns = c("CENT_X", "CENT_Y")) {
  data_sf <- st_as_sf(data, coords = coords_columns, crs = 4326)
  spatial_blocks <- cv_spatial(
    x = data_sf,
    k = cv,
    user_blocks = polygon,
    selection = "random",
    iteration = 100,
    progress = FALSE
  )
  spatial_blocks$folds_list
}

# --- 5. Optimization function ---
optimize_gbm <- function(response, lr, bf, tc, data, predictors, cv, spatial_folds) {
  formula <- as.formula(paste(response, "~", paste(predictors, collapse = "+")))
  cv_rmse <- numeric(cv)
  cv_trees <- numeric(cv)
  
  # Number of candidate trees: more for low learning rates
  n_candidates <- ifelse(lr <= 0.01, 10000, 5000)
  
  for (i in seq_len(cv)) {
    train_data <- data[spatial_folds[[i]][[1]], ]
    test_data  <- data[spatial_folds[[i]][[2]], ]
    
    if (nrow(train_data) == 0 || nrow(test_data) == 0) {
      stop("Empty fold detected: ", i)
    }
    
    model <- gbm(
      formula = formula, data = train_data,
      distribution = "gaussian", n.trees = n_candidates,
      shrinkage = lr, interaction.depth = tc, bag.fraction = bf,
      cv.folds = 0, n.minobsinnode = 10, verbose = FALSE
    )
    
    best_iter <- gbm.perf(model, plot.it = FALSE, method = "OOB")
    if (!is.numeric(best_iter) || length(best_iter) == 0) best_iter <- 100
    
    preds <- predict(model, test_data, n.tree = best_iter)
    cv_rmse[i]  <- sqrt(mean((preds - test_data[[response]])^2))
    cv_trees[i] <- best_iter
    gc()
  }
  
  data.frame(div = response, lr = lr, tc = tc, bf = bf,
             rmse = mean(cv_rmse), n_trees = round(mean(cv_trees)),
             stringsAsFactors = FALSE)
}

# --- 6. Run optimization for each response ---
set.seed(324)
optim_results <- list()

for (response in responses) {
  message("\nOptimizing for response: ", response)
  spatial_folds <- create_spatial_folds(df, polygon = polygon, cv = 10)
  
  response_result <- lapply(seq_len(nrow(hyper_grid)), function(i) {
    params <- hyper_grid[i, ]
    tryCatch({
      optimize_gbm(
        response = response, lr = params$lr, bf = params$bf, tc = params$tc,
        data = df, predictors = predictors, cv = 10,
        spatial_folds = spatial_folds
      )
    }, error = function(e) {
      message("Failed for ", response, " with params: ",
              paste(params, collapse = ", "), " - ", e$message)
      NULL
    })
  })
  optim_results[[response]] <- bind_rows(compact(response_result))
}

# --- 7. Select best parameters ---
final_results <- bind_rows(optim_results)
best_params <- final_results |>
  group_by(div) |>
  slice(which.min(rmse)) |>
  ungroup()

message("\nBest parameters:")
print(best_params)

# --- 8. Save ---
write.csv(best_params, file.path(output_dir, "brt_best_params.csv"), row.names = FALSE)
saveRDS(optim_results, file.path(output_dir, "brt_optimization_results.rds"))
message("03c complete.")
