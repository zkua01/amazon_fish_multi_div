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