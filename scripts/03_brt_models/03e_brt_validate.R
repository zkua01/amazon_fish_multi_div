# =============================================================================
# Script:  03e_brt_validate.R
# Purpose: Perform 10-fold spatial cross-validation with the optimized
#          hyperparameters. Computes out-of-sample R² (using training-fold
#          mean as reference), model RMSE, and null-model RMSE.
#
# Inputs:
#   data/processed/dataprepared_filtered.rds
#   output/intermediate/brt_best_params.csv
#
# Outputs:
#   output/intermediate/brt_cv_summary.csv
#
# Dependencies: gbm3, blockCV, sf, tidyverse, here
# =============================================================================

library(gbm3)
library(blockCV)
library(sf)
library(tidyverse)
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

# --- Helper functions ---
create_spatial_folds <- function(data, cv, polygon) {
  data_sf <- st_as_sf(data, coords = c("CENT_X", "CENT_Y"), crs = 4326)
  cv_spatial(x = data_sf, k = cv, user_blocks = polygon,
             selection = "random", iteration = 100, progress = FALSE)$folds_list
}

cv_results <- data.frame(
  response = character(), fold = integer(),
  model_rmse = numeric(), model_r2 = numeric(),
  null_rmse = numeric(), rmse_improvement_percent = numeric(),
  stringsAsFactors = FALSE
)

for (response in responses) {
  message("Validating response: ", response)
  params <- best_params |> filter(div == response)
  spatial_folds <- create_spatial_folds(df, cv = 10, polygon = polygon)
  
  for (i in 1:10) {
    train_data <- df[spatial_folds[[i]][[1]], ]
    test_data  <- df[spatial_folds[[i]][[2]], ]
    if (nrow(train_data) == 0 || nrow(test_data) == 0) next
    
    train_mean <- mean(train_data[[response]], na.rm = TRUE)
    
    model <- gbm(
      formula = as.formula(paste(response, "~", paste(predictors, collapse = "+"))),
      data = train_data, distribution = "gaussian",
      n.trees = params$n_trees, shrinkage = params$lr,
      interaction.depth = params$tc, bag.fraction = params$bf,
      cv.folds = 0, n.minobsinnode = 10, verbose = FALSE
    )
    
    preds  <- predict(model, test_data, n.trees = params$n_trees)
    actual <- test_data[[response]]
    
    model_rmse <- sqrt(mean((preds - actual)^2, na.rm = TRUE))
    ss_res <- sum((preds - actual)^2, na.rm = TRUE)
    ss_tot <- sum((actual - train_mean)^2, na.rm = TRUE)
    model_r2 <- 1 - (ss_res / ss_tot)
    
    null_rmse <- sqrt(mean((rep(train_mean, nrow(test_data)) - actual)^2, na.rm = TRUE))
    rmse_improvement <- ifelse(null_rmse > 0,
                               ((null_rmse - model_rmse) / null_rmse) * 100, 0)
    
    cv_results <- rbind(cv_results, data.frame(
      response = response, fold = i,
      model_rmse = model_rmse, model_r2 = model_r2,
      null_rmse = null_rmse, rmse_improvement_percent = rmse_improvement
    ))
  }
}

cv_summary <- cv_results |>
  group_by(response) |>
  summarise(across(c(model_rmse, model_r2, null_rmse, rmse_improvement_percent),
                   list(mean = mean, sd = sd), na.rm = TRUE),
            n_folds = max(fold), .groups = "drop")

write.csv(cv_summary, file.path(inter_dir, "brt_cv_summary.csv"), row.names = FALSE)
message("03e complete.")
