library(maps)
library(ggplot2)
library(RColorBrewer)
library(cowplot) 
library(gbm3)
library(Rmisc)
library(corrplot)
library(boot)
library(tidyverse)
library(purrr)
library(blockCV)
library(sf)
library(ggpubr)
library(ggtext)

path <- "Scripts/input/"
load(paste(path,'1.0b dataprepared_filt_fw',sep = ""))
basins_raw <- st_read("envdata/1.AmazonFish-distribution-subdrainage/SubDrainageShapefile/SubDrainageShapefile.shp")
basins <- basins_raw %>%
  st_transform(crs = 4326) %>%
  arrange(SBD_CD)

## BRT models: find optimal combination
dp<-na.omit(dataprepared_filt) %>% # remove NAs
  left_join(dplyr::select(basins, SBD_CD, geometry))
div<-dplyr::select(dp, muDiv, logSR, sesFD, sesPD) # select response variables
polygon<-dplyr::select(dp, SBD_CD, geometry) %>% st_as_sf() %>% st_make_valid() # select basin geometrics
coord <- dplyr::select(dp, CENT_X, CENT_Y) # select coordinates
env<-dplyr::select(dp, AREA:sample_area) %>% # select predictors
  mutate(across(c(AREA), log1p)) %>% #log+1 transformation on selected variables
  mutate_all(scale) %>% #standardize all variables apart from PebasLake (binary)
  mutate_all(as.numeric) #make sure variables are numeric and remove attributes

df <- bind_cols(div,env,coord)
namediv<-colnames(div)
predictors <- names(env)

hyper_grid <- expand.grid(
  lr = c(0.01,0.005,0.001),
  bf = seq(0.5,0.8,0.15),
  tc = seq(2,5,1)
  )

optim_results <- list()

# Create spatial blocks ONCE per response variable outside the optimization loop
create_spatial_folds <- function(data, cv, polygon, coords_columns = c("CENT_X", "CENT_Y")) {
  data_sf <- sf::st_as_sf(data, coords = coords_columns, crs = 4326)
  spatial_blocks <- blockCV::cv_spatial(
    x = data_sf,
    k = cv,
    user_blocks = polygon, 
    selection = "random",
    iteration = 100,
    progress = FALSE
  )
  return(spatial_blocks$folds_list)
}

optimize_gbm <- function(response, lr, bf, tc, data, predictors, cv, spatial_folds){
  formula = as.formula(paste(response, "~", paste(predictors, collapse = "+")))
  folds <- spatial_folds
  
  cv_rmse <- numeric(cv)
  cv_trees <- numeric(cv)
  
  for (i in 1:cv) {
    train_indices <- folds[[i]][[1]]
    test_indices <- folds[[i]][[2]]
    
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    # Check for valid data splits
    if(nrow(train_data) == 0 || nrow(test_data) == 0) {
      stop(paste("Empty fold detected:", i))
    }
    
    model <- gbm(
      formula = formula,
      data = train_data,
      distribution = "gaussian",
      n.trees = ifelse(lr<=0.01, 10000, 5000),
      shrinkage = lr,
      interaction.depth = tc,
      bag.fraction = bf,
      cv.folds = 0,
      n.minobsinnode = 10,
      verbose = F)
    
    best_iter <- gbm.perf(model, plot.it = F, method = "OOB")
    # Check for valid best iteration
    if(!is.numeric(best_iter) || length(best_iter) == 0) {
      message("No valid best iteration found, using default 100 trees")
      best_iter <- 100
    }
    
    preds <- predict(model, test_data, n.tree = best_iter)
    cv_rmse[i] <- sqrt(mean((preds - test_data[[response]])^2))
    cv_trees[i] <- best_iter
    
    gc()  # Memory cleanup
  }
  
  data.frame(
    div = response,
    lr = lr,
    tc = tc,
    bf = bf,
    rmse = mean(cv_rmse),
    #rmse_sd = sd(cv_rmse),
    n_trees = round(mean(cv_trees)),
    stringsAsFactors = F
  )
}

set.seed(324)
t0<-Sys.time()
for(response in namediv){
  message("\nOptimizing for response: ", response)
  spatial_folds <- create_spatial_folds(df, polygon = polygon, cv = 10)
  response_result <- lapply(1:nrow(hyper_grid), function(i){
    params <- hyper_grid[i,]
    tryCatch({
      optimize_gbm(
        response = response,
        lr = params$lr,
        bf = params$bf,
        tc = params$tc,
        data = df,
        predictors = predictors,
        cv = 10,
        spatial_folds = spatial_folds
      )
    }, error = function(e){
      message("Failed for ", response, " with params: ", paste(params, collapse = ", "), " - ", e$message)
      return(NULL)
    })
  })
  optim_results[[response]] <- bind_rows(compact(response_result))
}

final_results <- bind_rows(optim_results)
best_params <- final_results %>% 
  group_by(div) %>% 
  slice(which.min(rmse)) %>% ungroup()

# testing with cv=10 results in similar hyperparameters, with relatively small trees, suggesting that OOB is sufficient for optimization

#######################################

## BRT results

# Fit BRT model with optimized hyperparameters
fit_brt_model <- function(data, response, predictors, params){
  formula <- as.formula(paste(response,"~",paste(predictors, collapse = "+")))
  tryCatch({
    gbm(formula = formula,
        data = data,
        distribution = "gaussian",
        n.trees = params$n_trees,
        shrinkage = params$lr,
        interaction.depth = params$tc,
        bag.fraction = params$bf,
        cv.folds = 0,
        n.minobsinnode = 10,
        verbose = F)
  }, error = function(e){
    message("Model fit failed: ", e$message)
    return(NULL)
  })
}

# functions for bootstrapping

## relative influence
foo_influence <- function(data, indices, response, predictors, params){
  dt <- data[indices,]
  
  model <- fit_brt_model(dt, response, predictors, params)
  
  if(is.null(model)) return(rep(NA, length(predictors)))
  
  # get relative influence
  rel_inf <- relative_influence(model, num_trees = params$n_trees)
  (rel_inf[predictors] / sum(rel_inf[predictors])) * 100 # convert to percentage
}
## r2 values
foo_r2 <- function(data, indices, response, predictors, params){
  dt <- data[indices,]
  if(nrow(dt)<10) return (NA)
  model <- fit_brt_model(dt, response, predictors, params)
  if(is.null(model)) return(NA)
  
  # Get predictions on bootstrap sample
  preds <- predict(model, dt, n.trees = params$n_trees)
  actual <- dt[[response]]
  
  if(any(is.na(preds))||any(is.na(actual))) return(NA)
  
  # Calculate R²
  ss_res <- sum((preds - actual)^2, na.rm = TRUE)
  ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  if(ss_tot==0||is.na(ss_tot)) return(NA)
  r2 <- 1 - (ss_res / ss_tot)
  return(r2)
}
## pdps
foo_pdp <- function(data, indices, response, predictors, params, focal_var, grid_points = 50){
  dt <- data[indices,]
  
  model <- fit_brt_model(dt, response, predictors, params)
  if(is.null(model)) return(rep(NA, grid_points))
  
  # create grid for focal variable
  var_range <- range(dt[[focal_var]], na.rm=T)
  grid_vals <- seq(var_range[1], var_range[2], length.out=grid_points)
  
  # create reference data
  ref_data <- dt[1,] %>% slice(rep(1, length(grid_vals)))
  ref_data[[focal_var]] <- grid_vals
  
  # make predictions
  predict(model, newdata = ref_data, n.trees = params$n_trees)
}
## friedman's H
get_h_from_interact <- function(h_result) {
  if (is.null(h_result)) return(NA_real_)
  # atomic numeric scalar
  if (is.numeric(h_result) && length(h_result) == 1) return(as.numeric(h_result))
  # named numeric vector (e.g., c(H = 0.3) or other)
  if (is.numeric(h_result) && length(h_result) > 1) {
    if (!is.null(names(h_result)) && "H" %in% names(h_result)) return(as.numeric(h_result["H"]))
    return(as.numeric(h_result[1]))
  }
  # data.frame or tibble
  if (is.data.frame(h_result)) {
    if ("H" %in% colnames(h_result)) return(as.numeric(h_result$H[1]))
    # fallback: first numeric column/element
    numeric_cols <- sapply(h_result, is.numeric)
    if (any(numeric_cols)) return(as.numeric(h_result[[which(numeric_cols)[1]]][1]))
    return(NA_real_)
  }
  # list: look for H, else unlist
  if (is.list(h_result)) {
    if ("H" %in% names(h_result)) return(as.numeric(h_result$H))
    tmp <- unlist(h_result)
    if (length(tmp) > 0) return(as.numeric(tmp[1]))
    return(NA_real_)
  }
  # last resort
  return(NA_real_)
}

foo_friedman_h <- function(data, indices, response, predictors, params){
  dt <- data[indices, , drop=F]
  model <- fit_brt_model(dt,response, predictors, params)
  if(is.null(model)) return(rep(NA, choose(length(predictors), 2)))
  
  pairs <- combn(seq_along(predictors), 2, simplify = F)
  
  h_values <- numeric(length(pairs))
  
  for(j in seq_along(pairs)){
    pair <- pairs[[j]]
    tryCatch({
      h_result <- gbm3::interact(
        gbm_fit_obj = model,
        data = dt,
        var_indices = pair,
        num_trees = params$n_trees
      )
    }, error = function(e){
      message("H calculation failed for ", paste(pair, collapse = " x "), " : ", e$message)
      h_values[j] <- NA
    })
    h_values[j] <- get_h_from_interact(h_result)
  }
  return(h_values)
}

###################
# MANUAL SPATIAL BOOTSTRAP FUNCTION (simpler and more reliable)
manual_spatial_bootstrap <- function(data, response, predictors, params, 
                                     spatial_folds, stat_function, R = 1000, 
                                     additional_args = NULL) {
  
  n_folds <- length(spatial_folds)
  results <- list()
  
  cat("Running manual spatial bootstrap (", R, " iterations)...\n", sep = "")
  
  for (i in 1:R) {
    if (i %% 100 == 0) cat("  Iteration", i, "/", R, "\n")
    
    # Sample fold indices with replacement
    sampled_fold_indices <- sample(1:n_folds, size = n_folds, replace = TRUE)
    
    # Get training data from sampled folds
    train_indices <- unlist(lapply(sampled_fold_indices, function(idx) {
      spatial_folds[[idx]][[1]]  # Training indices for this fold
    }))
    
    # Create bootstrap sample
    bootstrap_data <- data[train_indices, ]
    
    # Skip if not enough data
    if (nrow(bootstrap_data) < 10) {
      next
    }
    
    # Calculate statistic with proper error handling
    tryCatch({
      if (!is.null(additional_args)) {
        result <- do.call(stat_function, c(list(
          data = bootstrap_data,
          response = response,
          predictors = predictors,
          params = params
        ), additional_args))
      } else {
        result <- stat_function(
          data = bootstrap_data,
          response = response,
          predictors = predictors,
          params = params
        )
      }
      results[[i]] <- result
    }, error = function(e) {
      # Skip failed iterations but continue
      message("  Iteration ", i, " failed: ", e$message)
    })
  }
  
  # Remove NULL results and convert to matrix/data.frame
  results <- results[!sapply(results, is.null)]
  
  if (length(results) == 0) {
    stop("All bootstrap iterations failed")
  }
  
  # Convert results to proper format
  if (is.list(results[[1]])) {
    # Handle list results (like PDPs)
    return(results)
  } else if (is.matrix(results[[1]]) || is.data.frame(results[[1]])) {
    # Handle matrix/data.frame results
    return(results)
  } else {
    # Handle vector results (like influence, R², H)
    result_matrix <- do.call(rbind, results)
    return(result_matrix)
  }
}

# REVISED STATISTIC FUNCTIONS (simplified signatures)
calc_influence <- function(data, response, predictors, params) {
  model <- fit_brt_model(data, response, predictors, params)
  if (is.null(model)) return(rep(NA, length(predictors)))
  
  rel_inf <- relative_influence(model, num_trees = params$n_trees)
  total_inf <- sum(rel_inf[predictors], na.rm = TRUE)
  
  if (total_inf == 0 || is.na(total_inf)) {
    return(rep(NA, length(predictors)))
  }
  
  (rel_inf[predictors] / total_inf) * 100
}

calc_r2 <- function(data, response, predictors, params) {
  if (nrow(data) < 10) return(NA)
  model <- fit_brt_model(data, response, predictors, params)
  if (is.null(model)) return(NA)
  
  preds <- predict(model, data, n.trees = params$n_trees)
  actual <- data[[response]]
  
  if (any(is.na(preds)) || any(is.na(actual))) return(NA)
  
  ss_res <- sum((preds - actual)^2, na.rm = TRUE)
  ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  
  if (ss_tot == 0 || is.na(ss_tot)) return(NA)
  
  1 - (ss_res / ss_tot)
}

calc_pdp <- function(data, response, predictors, params, focal_var, grid_points = 50) {
  model <- fit_brt_model(data, response, predictors, params)
  if (is.null(model)) return(rep(NA, grid_points))
  
  var_range <- range(data[[focal_var]], na.rm = TRUE)
  if (diff(var_range) == 0) {
    # Constant variable case
    ref_data <- data[1, ]
    pred <- predict(model, ref_data, n.trees = params$n_trees)
    return(rep(pred, grid_points))
  }
  
  grid_vals <- seq(var_range[1], var_range[2], length.out = grid_points)
  ref_data <- data[1, ] %>% slice(rep(1, length(grid_vals)))
  ref_data[[focal_var]] <- grid_vals
  
  predict(model, newdata = ref_data, n.trees = params$n_trees)
}

calc_friedman_h <- function(data, response, predictors, params) {
  model <- fit_brt_model(data, response, predictors, params)
  if (is.null(model)) return(rep(NA, choose(length(predictors), 2)))
  
  pairs <- combn(seq_along(predictors), 2, simplify = FALSE)
  h_values <- numeric(length(pairs))
  
  for (j in seq_along(pairs)) {
    pair <- pairs[[j]]
    tryCatch({
      h_result <- gbm3::interact(
        gbm_fit_obj = model,
        data = data,
        var_indices = pair,
        num_trees = params$n_trees
      )
      h_values[j] <- get_h_from_interact(h_result)
    }, error = function(e) {
      h_values[j] <- NA
    })
  }
  return(h_values)
}

# COMPLETELY REVISED RUN_ANALYSIS FUNCTION
run_analysis <- function(data, response, predictors, params_df, polygon,
                         n_influence_reps = 1000, n_pdp_reps = 200, n_h_reps = 200) {
  
  params <- params_df %>% filter(div == response)
  
  cat("\n=== Processing response variable:", response, "===\n")
  cat("  Using parameters: n_trees =", params$n_trees, 
      ", lr =", params$lr, 
      ", tc =", params$tc, 
      ", bf =", params$bf, "\n")
  
  # Create spatial folds ONCE
  spatial_folds <- create_spatial_folds(data, cv = 10, polygon = polygon)
  
  # ===== RELATIVE INFLUENCE WITH MANUAL SPATIAL BOOTSTRAP =====
  cat("\n1. Bootstrapping relative influence...\n")
  set.seed(324)
  influence_results <- manual_spatial_bootstrap(
    data = data,
    response = response,
    predictors = predictors,
    params = params,
    spatial_folds = spatial_folds,
    stat_function = calc_influence,
    R = n_influence_reps
  )
  
  # Calculate confidence intervals for influence
  influence_ci <- t(apply(influence_results, 2, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) < 10) return(c(NA, NA, NA))
    Rmisc::CI(x_clean)
  }))
  colnames(influence_ci) <- c("lower", "mean", "upper")
  rownames(influence_ci) <- predictors
  
  # ===== R² WITH MANUAL SPATIAL BOOTSTRAP =====
  cat("\n2. Bootstrapping R² value...\n")
  set.seed(324)
  r2_results <- manual_spatial_bootstrap(
    data = data,
    response = response,
    predictors = predictors,
    params = params,
    spatial_folds = spatial_folds,
    stat_function = calc_r2,
    R = n_influence_reps
  )
  
  # Calculate confidence interval for R²
  r2_clean <- r2_results[,1][!is.na(r2_results[,1])]
  if (length(r2_clean) < 10) {
    r2_ci <- c(NA, NA, NA)
  } else {
    r2_ci <- Rmisc::CI(r2_clean)
  }
  names(r2_ci) <- c("lower", "mean", "upper")
  
  # ===== PARTIAL DEPENDENCE PLOTS WITH MANUAL SPATIAL BOOTSTRAP =====
  cat("\n3. Calculating partial dependence plots...\n")
  
  pdp_list <- list()
  pdp_results_list <- list()
  
  for (var in predictors) {
    cat("  Variable:", var, "\n")
    
    set.seed(324)
    pdp_results <- manual_spatial_bootstrap(
      data = data,
      response = response,
      predictors = predictors,
      params = params,
      spatial_folds = spatial_folds,
      stat_function = calc_pdp,
      R = n_pdp_reps,
      additional_args = list(focal_var = var, grid_points = 50)
    )
    
    pdp_results_list[[var]] <- pdp_results
    
    # Create evaluation grid from ORIGINAL data
    var_range <- range(data[[var]], na.rm = TRUE)
    if (diff(var_range) == 0) {
      grid_vals <- rep(var_range[1], 50)
    } else {
      grid_vals <- seq(var_range[1], var_range[2], length.out = 50)
    }
    
    # Convert list of results to matrix
    pdp_matrix <- pdp_results
    
    # Calculate summary statistics with NA handling
    pdp_means <- apply(pdp_matrix, 2, function(x) mean(x, na.rm = TRUE))
    pdp_lower <- apply(pdp_matrix, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
    pdp_upper <- apply(pdp_matrix, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
    
    pdp_list[[var]] <- data.frame(
      variable = var,
      value = grid_vals,
      mean = pdp_means,
      lower = pdp_lower,
      upper = pdp_upper
    )
  }
  
  # ===== FRIEDMAN'S H WITH MANUAL SPATIAL BOOTSTRAP =====
  cat("\n4. Bootstrapping Friedman's H...\n")
  set.seed(324)
  h_results <- manual_spatial_bootstrap(
    data = data,
    response = response,
    predictors = predictors,
    params = params,
    spatial_folds = spatial_folds,
    stat_function = calc_friedman_h,
    R = n_h_reps
  )
  
  # Calculate confidence intervals for H statistics
  pairs <- combn(predictors, 2, simplify = FALSE)
  pair_names <- sapply(pairs, function(p) paste(p, collapse = " x "))
  
  h_ci <- t(apply(h_results, 2, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) < 10) return(c(NA, NA, NA))
    Rmisc::CI(x_clean)
  }))
  colnames(h_ci) <- c("lower", "mean", "upper")
  rownames(h_ci) <- pair_names
  
  # ===== RETURN RESULTS =====
  cat("\n✅ Successfully completed analysis for", response, "\n")
  
  list(
    response = response,
    parameters = params,
    influence = influence_ci,
    r2 = r2_ci,
    pdp = pdp_list,
    friedman_h = h_ci,
    bootstrap_raw = list(
      influence = influence_results,
      r2 = r2_results,
      pdp = pdp_results_list,
      friedman_h = h_results
    )
  )
}

# EXECUTION WITH PROPER ERROR HANDLING
t1 <- Sys.time()

response_var <- namediv
predictor_vars <- predictors

# Run analysis for each response
final_results <- list()
for (resp in response_var) {
  cat("\n\n", strrep("=", 50), "\n")
  cat("PROCESSING RESPONSE:", resp, "\n")
  cat(strrep("=", 50), "\n\n")
  
  tryCatch({
    result <- run_analysis(
      data = df,
      response = resp,
      predictors = predictor_vars,
      params_df = best_params,
      polygon = polygon,
      n_influence_reps = 1000,    # Start with fewer reps for testing
      n_pdp_reps = 200,           # Start with fewer reps for testing
      n_h_reps = 200              # Start with fewer reps for testing
    )
    final_results[[resp]] <- result
  }, error = function(e) {
    cat("\n❌ ERROR processing", resp, ":", e$message, "\n")
    cat("Error occurred at:", Sys.time(), "\n")
    final_results[[resp]] <- NULL
  })
}

final_results <- final_results[!sapply(final_results, is.null)]

t2 <- Sys.time()
cat("\n\n", strrep("=", 50), "\n")
cat("TOTAL PROCESSING TIME:", t2 - t1, "\n")
cat(strrep("=", 50), "\n")
###################
saveRDS(final_results, paste0(path,"brt_final_results_fw.rds"))

final_results <- read_rds(paste0(path,"brt_final_results_fw.rds"))

# model validation
# switch responses for all 4 validations
cv_results <- data.frame(
  response = character(),
  fold = integer(),
  model_rmse = numeric(),
  model_r2 = numeric(),
  null_rmse = numeric(),
  rmse_improvement_percent = numeric(),
  stringsAsFactors = FALSE
)
for (response in namediv) {
  cat("\nProcessing response variable:", response, "\n")
  
  # Get best parameters for this response
  params <- final_results[[response]]$parameters
  if (nrow(params) == 0) {
    stop(paste("No best parameters found for response:", response))
  }
  
  # Create spatial folds for this response
  spatial_folds <- create_spatial_folds(df, cv = 10, polygon = polygon)
  
  # Process each fold
  for (i in 1:10) {
    cat("  Fold", i, "of 10...\n")
    
    train_indices <- spatial_folds[[i]][[1]]
    test_indices <- spatial_folds[[i]][[2]]
    
    train_data <- df[train_indices, ]
    test_data <- df[test_indices, ]
    
    # Check for valid data splits
    if (nrow(train_data) == 0 || nrow(test_data) == 0) {
      warning(paste("Empty fold detected for response", response, "fold", i))
      next
    }
    
    # Get training mean for proper CV R² calculation
    train_mean <- mean(train_data[[response]], na.rm = TRUE)
    
    # Train model with optimal parameters
    model <- gbm(
      formula = as.formula(paste(response, "~", paste(predictors, collapse = "+"))),
      data = train_data,
      distribution = "gaussian",
      n.trees = params$n_trees,
      shrinkage = params$lr,
      interaction.depth = params$tc,
      bag.fraction = params$bf,
      cv.folds = 0,
      n.minobsinnode = 10,
      verbose = FALSE
    )
    
    # Predict on test fold
    preds <- predict(model, test_data, n.trees = params$n_trees)
    actual <- test_data[[response]]
    
    # Calculate model metrics
    model_rmse <- sqrt(mean((preds - actual)^2, na.rm = TRUE))
    
    # Calculate cross-validated R² (using training mean for proper out-of-sample calculation)
    ss_res <- sum((preds - actual)^2, na.rm = TRUE)
    ss_tot <- sum((actual - train_mean)^2, na.rm = TRUE)
    model_r2 <- 1 - (ss_res / ss_tot)
    
    # Calculate null model metrics (using training mean as prediction)
    null_preds <- rep(train_mean, nrow(test_data))
    null_rmse <- sqrt(mean((null_preds - actual)^2, na.rm = TRUE))
    
    # Calculate RMSE improvement percentage
    rmse_improvement_percent <- ifelse(null_rmse > 0, 
                                       ((null_rmse - model_rmse) / null_rmse) * 100, 
                                       0)
    
    # Store results
    fold_result <- data.frame(
      response = response,
      fold = i,
      model_rmse = model_rmse,
      model_r2 = model_r2,
      null_rmse = null_rmse,
      rmse_improvement_percent = rmse_improvement_percent
    )
    
    cv_results <- rbind(cv_results, fold_result)
    
    # Memory cleanup
    rm(model, preds, actual, null_preds)
    gc()
  }
}

# Calculate summary statistics by response
cv_summary <- cv_results %>%
  group_by(response) %>%
  dplyr::summarise(model_rmse_mean = mean(model_rmse, na.rm = TRUE),
            model_rmse_sd = sd(model_rmse, na.rm = TRUE),
            model_r2_mean = mean(model_r2, na.rm = TRUE),
            model_r2_sd = sd(model_r2, na.rm = TRUE),
            null_rmse_mean = mean(null_rmse, na.rm = TRUE),
            null_rmse_sd = sd(null_rmse, na.rm = TRUE),
            rmse_improvement_mean = mean(rmse_improvement_percent, na.rm = TRUE),
            rmse_improvement_sd = sd(rmse_improvement_percent, na.rm = TRUE),
            n_folds = max(fold), .groups = "drop")

as.data.frame(cv_summary)
saveRDS(cv_summary, paste0(path,"brt_cv_summary_fw.rds"))

# save final results influence tables
final_results$logSR$influence |>
  as.data.frame() |>
  rownames_to_column(var = "predictor") |>
  arrange(-mean) |>
  write.csv(paste0(path,"brt_influence_td_fw.csv"), row.names = FALSE)
final_results$sesPD$influence |>
  as.data.frame() |>
  rownames_to_column(var = "predictor") |>
  arrange(-mean) |>
  write.csv(paste0(path,"brt_influence_pd_fw.csv"), row.names = FALSE)
final_results$sesFD$influence |>
  as.data.frame() |>
  rownames_to_column(var = "predictor") |>
  arrange(-mean) |>
  write.csv(paste0(path,"brt_influence_fd_fw.csv"), row.names = FALSE)
final_results$muDiv$influence |>
  as.data.frame() |>
  rownames_to_column(var = "predictor") |>
  arrange(-mean) |>
  write.csv(paste0(path,"brt_influence_md_fw.csv"), row.names = FALSE)

#save parameters
bind_rows(final_results$logSR$parameters, final_results$sesPD$parameters, final_results$sesFD$parameters, final_results$muDiv$parameters) |>
  write.csv(paste0(path, "brt_parameters.csv"), row.names = FALSE)

# plots
## create categories for data frame
name<-data.frame(col=c(rep("Historical variables",3),
                       rep("Environmental variables",3),
                       rep('Habitat size',4),
                       rep('Positioning',1),
                       rep('Fragmentation',2),
                       rep('Effort', 1)),
                 abb=c('PC1_CurrPast', 'slr25m', 'slr100m',
                       'PC2_Temp', 'PC1_Water', 'PC4_Energy',
                       'AREA', 'NetworkDensity', 'soilDiv', 'LCDiv',
                       'DistMouth_km',
                       'PC1_Elev', 'waterfalls',
                       'sample_area'),
                 label = c('Hist. climate warming', 'SLR 25m', 'SLR 100m',
                           'Solar energy input', 'Water scarcity', 'Productivity seasonality',
                           'Area', 'Network density', 'Soil diversity', 'Land cover diversity',
                           'Distance to mouth (km)',
                           'Terrain Elevation', 'Waterfall counts',
                           'Sampling effort'))
lev<-c('Historical variables', 'Environmental variables', 'Anthropogenic variables', 'Habitat size', 'Positioning', 'Fragmentation', 'Effort')
name$col<-factor(name$col,levels =lev )
# set color palette
col <- c('#D55E00', '#117733', '#FF9999', '#9999FF', '#66B2FF', 'gray40')

plot_influence <- function(results_list, response_name){
  influence_data <- as.data.frame(results_list[[response_name]]$influence) %>%
    rownames_to_column("predictor") %>%
    left_join(name, by = c("predictor" = "abb")) %>%
    group_by(col) %>%
    dplyr::arrange(mean, .by_group = T) %>%
    ungroup()
  
  response_label = case_when(response_name == "muDiv" ~ "μ-diversity",
                             response_name == "sesPD" ~ "Phylogenetic diversity",
                             response_name == "sesFD" ~ "Functional diversity",
                             response_name == "logSR" ~ "Taxonomic diversity")
  
  params <- results_list[[response_name]]$parameters
  r2 <- results_list[[response_name]]$r2
  
  ggplot(influence_data, aes(x = reorder(label, mean), y = mean))+
    geom_linerange(aes(ymin = lower, ymax = upper, colour = col), size = 4, alpha = 0.5)+
    geom_point(aes(colour = col), size = 2.5)+
    geom_hline(yintercept = 7.14, color = "red", linetype = 2)+
    scale_colour_manual(name = "", values = col)+
    scale_fill_manual(name = "", values = col)+
    coord_flip(ylim = c(0,60))+
    labs(title = response_label,
         # subtitle = sprintf("%d trees, lr=%.3f, tc=%d, bf=%.2f, R<sup>2</sup>=%.2f", 
         #                    params$n_trees, params$lr, params$tc, params$bf, r2['mean']),
         x = NULL,
         y = "Relative influence (%)")+
    theme_classic()+
    theme(axis.text = ggplot2::element_text(size = 12, color = "black"),
          panel.grid = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.ticks.length.x = ggplot2::unit(.15, "cm"),
          legend.position = c(.9, 0.6),
          legend.background = ggplot2::element_rect(fill = "transparent"),
          panel.border = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_line(size = 0.5, colour = "black"),
          plot.margin = ggplot2::unit(c(2, 5, 0, 0), "mm"),
          plot.title = ggplot2::element_text(hjust = 0.5),
          axis.title = ggplot2::element_text(size = 16),
          plot.subtitle = ggtext::element_markdown(size = 8))
}

f1 <- plot_influence(results_list = final_results, response_name = response_var[1])
f2 <- plot_influence(results_list = final_results, response_name = response_var[2]) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
f3 <- plot_influence(results_list = final_results, response_name = response_var[3])
f4 <- plot_influence(results_list = final_results, response_name = response_var[4]) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

##bar plot 
influence_dfs = list()
for (i in seq_along(response_var)) {
  var_name <- final_results[[i]]$response
  influence_dfs[[var_name]] <- as.data.frame(final_results[[i]]$influence) %>%
    rownames_to_column("predictor") %>%
    left_join(name, by=c("predictor"="abb")) %>%
    dplyr::rename(var = col) %>%
    dplyr::group_by(var) %>%
    dplyr::summarise(percent = sum(mean)) %>%
    ungroup()
}

bar <- bind_rows(influence_dfs, .id = "div") |>
  dplyr::mutate(div = factor(div, levels = c("muDiv","logSR","sesPD","sesFD")),
                var = factor(var, levels = lev))

bplot<-function(x,y="none"){
  bar_sub<-bar[bar$div==x,]
  ggplot(bar_sub,aes(x=div,y=percent,fill=var))+
    geom_bar(stat="identity")+
    guides(fill = guide_legend(reverse=F))+
    scale_y_continuous(expand=c(0,0),limits = c(0,100.1), breaks = seq(0,100,25),labels = seq(0,100,25))+
    scale_fill_manual(values=col)+
    theme(panel.background = ggplot2::element_rect(fill = "transparent"),
          plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
          panel.grid = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(),
          axis.line.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(size = 10, color = "black"),
          axis.title = ggplot2::element_blank(), # Combines xlab("") and ylab("")
          legend.position = y,
          legend.title = ggplot2::element_blank(),
          legend.background = ggplot2::element_rect(fill = "transparent"),
          legend.box.background = ggplot2::element_rect(fill = "transparent", color = "transparent"),
          plot.margin = ggplot2::margin(5, 2, 2, 2, "mm"),
          plot.title = ggplot2::element_text(size = 18, hjust = 0.5))
}
b1<-bplot("muDiv")
b2<-bplot("logSR")
b4<-bplot("sesPD")
b3<-bplot("sesFD")

# get y-axis range from f plots
y1_range <- layer_scales(f1)$y$range$range
y2_range <- layer_scales(f2)$y$range$range
y3_range <- layer_scales(f3)$y$range$range
y4_range <- layer_scales(f4)$y$range$range

c1<-f1+annotation_custom(ggplotGrob(b1),ymin = 60*0.7,ymax = 60*0.9,xmin = 1, xmax = 10)
c2<-f2+annotation_custom(ggplotGrob(b2),ymin = 60*0.7,ymax = 60*0.9,xmin = 1, xmax = 10)
c3<-f3+annotation_custom(ggplotGrob(b3),ymin = 60*0.7,ymax = 60*0.9,xmin = 1, xmax = 10)
c4<-f4+annotation_custom(ggplotGrob(b4),ymin = 60*0.7,ymax = 60*0.9,xmin = 1, xmax = 10)

fig<-ggpubr::ggarrange(c2,c4,c3,c1, nrow = 2, ncol = 2, common.legend = T, legend = "bottom", labels = "AUTO")
print(fig)

ggsave(fig, filename = paste0("Scripts/input/Figures/BRT_fw_", format(Sys.Date(), "%Y%m%d"), ".png"), 
       width=5000,height=4000, units = "px", dpi = 400)

dev.off()

################################################
# plot partial dependency plots
## prep scale info for back transformation
env_trans<-dplyr::select(dp, AREA:sample_area) %>% # select predictors
  mutate(across(c(AREA), log1p)) %>% #log+1 transformation on selected variables
  mutate_all(scale)

plot_response_pdps <- function(results_list, response_name, ncol = 4, ref_data, log1p_vars = "AREA", use_original_scale=TRUE){
  
  #combine all PDP data for this response
  pdp_data <- bind_rows(results_list[[response_name]]$pdp, .id = "predictor") %>%
    mutate(predictor = factor(predictor, levels = names(results_list[[response_name]]$pdp)))
  
  #Get influence values for ordering
  influence_vals <- as.data.frame(results_list[[response_name]]$influence) %>%
    rownames_to_column("predictor") %>%
    dplyr::select(predictor, mean_influence = mean)
  
  #Extract transformation parameters from refdata
  pred_vars <- names(results_list[[response_name]]$pdp)
  trans_params <- tibble(predictor = pred_vars) %>%
    mutate(center = map_dbl(predictor, ~{a <- attr(ref_data[[.x]], "scaled:center"); if(is.null(a)) 0 else a}),
           scale = map_dbl(predictor, ~{a <- attr(ref_data[[.x]], "scaled:scale"); if(is.null(a)) 1 else a}),
           is_scaled = map_lgl(predictor, ~!is.null(attr(ref_data[[.x]], "scaled:center"))),
           is_log1p = predictor %in% log1p_vars)
  
  #merge influence value, order predictors, and JOIN LABELS
  pdp_data_plot <- pdp_data %>%
    left_join(influence_vals, by = "predictor") %>%
    left_join(trans_params, by = "predictor") %>%
    left_join(name %>% dplyr::select(abb, label), by = c("predictor" = "abb")) %>% # Added join
    mutate(predictor = fct_reorder(predictor, -mean_influence),
           label = fct_reorder(label, -mean_influence), # Reorder labels to match influence
           value_unscaled = ifelse(is_scaled, value * scale + center, value),
           value_orig = ifelse(is_log1p, expm1(value_unscaled), value_unscaled)) %>%
    dplyr::select(-center,-scale,-is_scaled,-is_log1p,-value_unscaled)
  
  response_label = case_when(response_name == "muDiv" ~ "μ-diversity",
                             response_name == "sesPD" ~ "Phylogenetic diversity",
                             response_name == "sesFD" ~ "Functional diversity",
                             response_name == "logSR" ~ "Taxonomic diversity")
  #x-axis column and label
  x_col <- if(use_original_scale) "value_orig" else "value"
  
  #create plot using `label` for faceting
  ggplot(pdp_data_plot, aes(x = .data[[x_col]], y = mean))+
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'darkorange', alpha = 0.2)+
    geom_line(color = 'darkorange', linewidth = 0.8)+
    facet_wrap(~label, scales = "free_x", ncol = ncol)+ # Changed from ~predictor to ~label
    labs(title = paste("Partial dependence plots: ", response_label),
         # subtitle = paste("n_trees = ", results_list[[response_name]]$parameters$n_trees,
         #                  "lr = ", results_list[[response_name]]$parameters$lr,
         #                  "tc = ", results_list[[response_name]]$parameters$tc,
         #                  "bf = ", results_list[[response_name]]$parameters$bf),
         x = "Predictor value",
         y = "Partial dependence")+
    theme_classic()+
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_rect(fill = "grey90", color = NA),
          panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 9, color = "grey30"))
}

plot_response_pdps_t5 <- function(results_list, response_name, ref_data, log1p_vars = "AREA", ncol = 4, use_original_scale = TRUE, return_data = FALSE){
  
  #combine all PDP data for this response
  pdp_data <- bind_rows(results_list[[response_name]]$pdp, .id = "predictor") %>%
    mutate(predictor = factor(predictor, levels = names(results_list[[response_name]]$pdp)))
  
  #Get influence values for ordering
  influence_vals <- as.data.frame(results_list[[response_name]]$influence) %>%
    rownames_to_column("predictor") %>%
    dplyr::select(predictor, mean_influence = mean)
  
  #extract transformation parameters
  pred_vars <- names(results_list[[response_name]]$pdp)
  trans_params <- tibble(predictor = pred_vars) %>%
    mutate(center = map_dbl(predictor, ~{a <- attr(ref_data[[.x]], "scaled:center"); if(is.null(a)) 0 else a}),
           scale = map_dbl(predictor, ~{a <- attr(ref_data[[.x]], "scaled:scale"); if(is.null(a)) 1 else a}),
           is_scaled = map_lgl(predictor, ~!is.null(attr(ref_data[[.x]], "scaled:center"))),
           is_log1p = predictor %in% log1p_vars)
  
  influence_t5 <- head(arrange(influence_vals, desc(mean_influence)), 6)$predictor
  
  #merge influence value, order predictors, and JOIN LABELS
  pdp_data_plot <- pdp_data %>%
    left_join(influence_vals, by = "predictor") %>%
    left_join(trans_params, by = "predictor") %>%
    left_join(name %>% dplyr::select(abb, label), by = c("predictor" = "abb")) %>% # Added join
    mutate(predictor = fct_reorder(predictor, -mean_influence),
           label = fct_reorder(label, -mean_influence)) %>%
    arrange(desc(mean_influence)) %>%
    filter(predictor %in% influence_t5) %>% # Fixed bug: changed 'variable' to 'predictor'
    mutate(value_unscaled = ifelse(is_scaled, value*scale+center, value),
           value_orig = ifelse(is_log1p, expm1(value_unscaled), value_unscaled)) %>%
    dplyr::select(-center,-scale,-is_scaled,-is_log1p,-value_unscaled)
  
  response_label = case_when(response_name == "muDiv" ~ "μ-diversity",
                             response_name == "sesPD" ~ "Phylogenetic diversity",
                             response_name == "sesFD" ~ "Functional diversity",
                             response_name == "logSR" ~ "Taxonomic diversity")
  
  x_col = if(use_original_scale) "value_orig" else "value"
  
  #create plot using `label` for faceting
  p <- ggplot(pdp_data_plot, aes(x = .data[[x_col]], y = mean))+
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'darkorange', alpha = 0.2)+
    geom_line(color = 'darkorange', linewidth = 0.8)+
    facet_wrap(~label, scales = "free_x", ncol = ncol)+ # Changed from ~predictor to ~label
    labs(title = response_label,
         # subtitle = paste("n_trees = ", results_list[[response_name]]$parameters$n_trees,
         #                  "lr = ", results_list[[response_name]]$parameters$lr,
         #                  "tc = ", results_list[[response_name]]$parameters$tc,
         #                  "bf = ", results_list[[response_name]]$parameters$bf),
         x = "Predictor value",
         y = "Partial dependence")+
    theme_classic()+
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_rect(fill = "grey90", color = NA),
          panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 9, color = "grey30"))
  
  if (return_data){
    return(list(plot = p, data = pdp_data_plot))
  } else {
    return(p)
  }
}

walk(names(final_results), function(resp){
  p <- plot_response_pdps(final_results, resp, ref_data = env_trans, ncol = 4)
  
  ggsave(filename = paste0(resp,"_pdps_fw.png"), plot = p, width = 5000, height = 4000, units = "px", dpi = 400)
  
  print(p)
})

md_pdp_t5 <- plot_response_pdps_t5(final_results, "muDiv", ref_data = env_trans, ncol = 3, return_data = TRUE)
td_pdp_t5 <- plot_response_pdps_t5(final_results, "logSR", ref_data = env_trans, ncol = 3, return_data = TRUE)
fd_pdp_t5 <- plot_response_pdps_t5(final_results, "sesFD", ref_data = env_trans, ncol = 3, return_data = TRUE)
pd_pdp_t5 <- plot_response_pdps_t5(final_results, "sesPD", ref_data = env_trans, ncol = 3, return_data = TRUE)

pdp_t5_plot <- ggarrange(td_pdp_t5$plot, pd_pdp_t5$plot, fd_pdp_t5$plot, md_pdp_t5$plot, ncol = 2, nrow = 2, labels = "AUTO")

pdp_t5_plot

ggsave(paste0("Scripts/input/Figures/BRT_PDP_T5_fw_", format(Sys.Date(), "%Y%m%d"), ".png"), pdp_t5_plot,
       width = 5000, height = 4000, units = "px", dpi = 400)

## Friedman's H
plot_friedman_h <- function(results_list, response_name){
  response_label = case_when(response_name == "muDiv" ~ "μ-diversity",
                             response_name == "sesPD" ~ "Phylogenetic diversity",
                             response_name == "sesFD" ~ "Functional diversity",
                             response_name == "logSR" ~ "Taxonomic diversity")
  
  # Create a named vector to map abbreviations to labels for string replacement
  replace_dict <- setNames(name$label, name$abb)
  
  #combine all H data for this response
  h_data <- as_tibble(results_list[[response_name]]$friedman_h, rownames = "pair") %>%
    arrange(desc(mean)) %>%
    mutate(pair = str_replace_all(pair, " x ", " : ")) %>%
    mutate(pair = str_replace_all(pair, replace_dict)) # Replaces abbreviations with labels
  
  h_data_t5 <- head(h_data, 6)
  
  #create plot
  ggplot(h_data_t5, aes(x = pair, y = mean))+
    geom_linerange(aes(ymin = lower, ymax = upper), color = 'navy', alpha = 0.2, size = 4)+
    geom_point(color = 'navy', size = 3)+
    geom_hline(aes(yintercept = 0.2), linetype = 2, colour = 'red')+
    labs(title = response_label,
         x = "Predictor pairs",
         y = "Friedman's H")+
    theme_classic()+
    coord_flip(ylim = c(0,0.5))+
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_rect(fill = "grey90", color = NA),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size=10,color = "black"),
          axis.text.x = element_text(size=10, colour = "black"),
          axis.title = element_text(size=14),
          title = element_text(size = 15))
}
# walk(names(final_results), function(resp){
#   p <- plot_friedman_h(final_results, resp)
#   
#   ggsave(filename = paste0(resp,"_freidman_h.png"), plot = p, width = 5000, height = 4000, units = "px", dpi = 400)
#   
#   print(p)
# })

plot_list <- purrr::map(names(final_results), function(resp) {
  plot_friedman_h(final_results, resp)
})

# Use ggarrange to combine all plots into one arranged figure
combined_hplot <- ggarrange(plotlist = list(plot_list[[2]], plot_list[[4]], plot_list[[3]], plot_list[[1]]), ncol = 2, nrow = 2, align = "hv") %>%
  annotate_figure(top = text_grob("Top 6 pairwise Friedman's H plots", size = 16))
combined_hplot

ggsave("Scripts/input/Figures/brt_t6_friedman_h_fw_plot.png", combined_hplot, width = 5000, height = 4000, units = 'px', dpi = 400)

library(gbm3)
library(iml)
library(ggpubr)

# 1. Create a dictionary to map abbreviations to your descriptive labels
label_dict <- setNames(name$label, name$abb)
label_dict_title <- label_dict
label_dict_title["DistMouth_km"] <- "Dist. mouth"
label_dict_title["PC1_CurrPast"] <- "Hist. climate"
label_dict_title["PC4_Energy"] <- "Prod. seasonality"
label_dict_title["PC2_Temp"] <- "Solar input"

# 2. Extract Friedman's H and split the pair strings into var1 and var2
h_table <- purrr::map_dfr(names(final_results), function(resp){
  h_mat <- final_results[[resp]]$friedman_h
  h_df <- as.data.frame(h_mat) %>%
    rownames_to_column("pair") %>%
    mutate(response = resp)
  return(h_df)
})

h_table_split <- h_table %>%
  tidyr::separate(pair, into = c("var1", "var2"), sep = " x ") %>%
  filter(mean > 0.2) # Criterion 1: Friedman's H > 0.2

# 3. Extract relative influence for all individual variables
inf_table <- purrr::map_dfr(names(final_results), function(resp) {
  inf_df <- as.data.frame(final_results[[resp]]$influence) %>%
    rownames_to_column("predictor") %>%
    mutate(response = resp)
  return(inf_df)
})

# 4. Filter pairs where BOTH variables have relative influence > 10%
# (Note: If you prefer the SUM of their influences to be > 10%, change the filter to: filter(inf1 + inf2 > 10))
valid_pairs <- h_table_split %>%
  left_join(inf_table %>% dplyr::select(response, predictor, inf1 = mean), by = c("response", "var1" = "predictor")) %>%
  left_join(inf_table %>% dplyr::select(response, predictor, inf2 = mean), by = c("response", "var2" = "predictor")) %>%
  filter(inf1 > 7.14 | inf2 > 7.14) %>% # Criterion 2: Relative influence > 10% for both
  dplyr::select(response, var1, var2, mean)

# 5. Refit models and create iml Predictor objects dynamically
models_list <- list()
preds_list <- list()

for (resp in names(final_results)) {
  params <- final_results[[resp]]$parameters
  
  # Refit model
  model_obj <- gbm(
    formula = as.formula(paste(resp, "~", paste(predictors, collapse = "+"))),
    data = df,
    distribution = "gaussian",
    n.trees = params$n_trees,
    shrinkage = params$lr,
    interaction.depth = params$tc,
    bag.fraction = params$bf,
    cv.folds = 0,
    n.minobsinnode = 10,
    verbose = F,
    train.fraction = 1.0,
    keep.data = T
  )
  models_list[[resp]] <- model_obj
  
  # Create iml Predictor
  preds_list[[resp]] <- iml::Predictor$new(
    model = model_obj,
    data = df,
    y = df[[resp]],
    predict.function = function(object, newdata) predict(object, newdata = newdata, n.trees = params$n_trees)
  )
}

# 6. Generate ICE and 2D PDP plots dynamically
plot_list_ice <- list()
plot_list_2dpdp <- list()

# Helper for response labels
get_resp_label <- function(resp) {
  case_when(
    resp == "muDiv" ~ "μ-diversity",
    resp == "sesPD" ~ "Phylogenetic diversity",
    resp == "sesFD" ~ "Functional diversity",
    resp == "logSR" ~ "Taxonomic diversity"
  )
}

# Base theme used in your script
base_theme <- theme_classic() +
  theme(panel.grid = element_blank(), axis.line = element_line(),
        panel.background = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 16), plot.margin = margin(5, 2, 2, 2, "mm"),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.text = element_text(size = 16))

for (resp in names(final_results)) {
  pairs_resp <- valid_pairs %>% filter(response == resp)
  
  if (nrow(pairs_resp) == 0) {
    message(paste("No valid pairs found for", resp))
    next
  }
  
  pred_obj <- preds_list[[resp]]
  resp_label <- get_resp_label(resp)
  
  # --- ICE Plots ---
  vars_ice <- unique(c(pairs_resp$var1, pairs_resp$var2))
  
  ice_obj <- iml::FeatureEffects$new(
    pred_obj,
    feature = vars_ice,
    method = "ice"
  )
  
  # Extract and center data dynamically using .id to track the feature
  ice_data <- bind_rows(ice_obj$results, .id = ".feature") %>%
    group_by(.feature, .id) %>%
    mutate(centered = .value - first(.value)) %>%
    ungroup()
  
  p_ice <- ggplot(ice_data, aes(x = .borders, y = centered, group = .id)) +
    geom_line(alpha = .3) +
    stat_summary(aes(group = 1), fun = mean, geom = "line", color = "red", size = 1) +
    labs(y = paste("Centered predicted", resp_label), x = "Predictor") +
    # Use labeller to automatically swap abbreviations for descriptive labels in the facets
    facet_wrap(~.feature, scales = "free", labeller = labeller(.feature = label_dict)) +
    base_theme
  
  plot_list_ice[[resp]] <- p_ice
  
  # --- 2D PDP Plots ---
  pdp_plots <- list()
  for (i in 1:nrow(pairs_resp)) {
    v1 <- pairs_resp$var1[i]
    v2 <- pairs_resp$var2[i]
    
    pdp_obj <- iml::FeatureEffect$new(
      pred_obj,
      feature = c(v1, v2),
      method = "pdp"
    )
    
    # Override axes and title with descriptive labels
    p_pdp <- plot(pdp_obj) +
      labs(x = label_dict[v1], y = label_dict[v2], 
           title = paste(label_dict_title[v1], "&", label_dict_title[v2]),
           fill = "Partial dependence") +
      base_theme
    
    pdp_plots[[i]] <- p_pdp
  }
  
  # Combine 2D PDPs for this response
  p_combined <- ggarrange(plotlist = pdp_plots, common.legend = TRUE, legend = "right", labels = "AUTO")
  plot_list_2dpdp[[resp]] <- p_combined
}

# 7. Save the plots
resp_names <- c("muDiv" = "mudiv", "logSR" = "td", "sesPD" = "pd", "sesFD" = "fd")

for (resp in names(final_results)) {
  if (!is.null(plot_list_ice[[resp]])) {
    ggsave(paste0(path, "Figures/", resp_names[resp], "_fw_ice_plots.png"), 
           plot_list_ice[[resp]], width = 5000, height = 4000, units = "px", dpi = 400)
  }
  if (!is.null(plot_list_2dpdp[[resp]])) {
    ggsave(paste0(path, "Figures/", resp_names[resp], "_fw_2dpdp_plot.png"), 
           plot_list_2dpdp[[resp]], width = 5000, height = 4000, units = "px", dpi = 400)
  }
}
