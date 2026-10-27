# =============================================================================
# Shared BRT utility functions for scripts 03d–03f
# =============================================================================

# Fit a BRT model with given parameters (returns NULL on failure)
fit_brt_model <- function(data, response, predictors, params) {
  formula <- as.formula(paste(response, "~", paste(predictors, collapse = "+")))
  tryCatch({
    gbm(formula = formula, data = data,
        distribution = "gaussian", n.trees = params$n_trees,
        shrinkage = params$lr, interaction.depth = params$tc,
        bag.fraction = params$bf, cv.folds = 0,
        n.minobsinnode = 10, verbose = FALSE)
  }, error = function(e) {
    message("Model fit failed: ", e$message)
    NULL
  })
}

# Extract scalar Friedman's H from various gbm3::interact return types
get_h_from_interact <- function(h_result) {
  if (is.null(h_result)) return(NA_real_)
  if (is.numeric(h_result) && length(h_result) == 1) return(as.numeric(h_result))
  if (is.numeric(h_result) && length(h_result) > 1) {
    if (!is.null(names(h_result)) && "H" %in% names(h_result)) {
      return(as.numeric(h_result["H"]))
    }
    return(as.numeric(h_result[1]))
  }
  if (is.data.frame(h_result)) {
    if ("H" %in% colnames(h_result)) return(as.numeric(h_result$H[1]))
    numeric_cols <- vapply(h_result, is.numeric, logical(1))
    if (any(numeric_cols)) return(as.numeric(h_result[[which(numeric_cols)[1]]][1]))
    return(NA_real_)
  }
  if (is.list(h_result)) {
    if ("H" %in% names(h_result)) return(as.numeric(h_result$H))
    tmp <- unlist(h_result)
    if (length(tmp) > 0) return(as.numeric(tmp[1]))
    return(NA_real_)
  }
  NA_real_
}

# --- Bootstrap statistic functions ---

calc_influence <- function(data, response, predictors, params) {
  model <- fit_brt_model(data, response, predictors, params)
  if (is.null(model)) return(rep(NA, length(predictors)))
  # Requires dismo::relative_influence or gbm3 equivalent
  rel_inf <- dismo::relative_influence(model, num_trees = params$n_trees)
  total_inf <- sum(rel_inf[predictors], na.rm = TRUE)
  if (total_inf == 0 || is.na(total_inf)) return(rep(NA, length(predictors)))
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
    ref_data <- data[1, ]
    return(rep(predict(model, ref_data, n.trees = params$n_trees), grid_points))
  }
  grid_vals <- seq(var_range[1], var_range[2], length.out = grid_points)
  ref_data <- data[1, ] |> slice(rep(1, length(grid_vals)))
  ref_data[[focal_var]] <- grid_vals
  predict(model, newdata = ref_data, n.trees = params$n_trees)
}

calc_friedman_h <- function(data, response, predictors, params) {
  model <- fit_brt_model(data, response, predictors, params)
  if (is.null(model)) return(rep(NA, choose(length(predictors), 2)))
  pairs <- combn(seq_along(predictors), 2, simplify = FALSE)
  h_values <- numeric(length(pairs))
  for (j in seq_along(pairs)) {
    tryCatch({
      h_result <- gbm3::interact(
        gbm_fit_obj = model, data = data,
        var_indices = pairs[[j]], num_trees = params$n_trees
      )
      h_values[j] <- get_h_from_interact(h_result)
    }, error = function(e) { h_values[j] <- NA })
  }
  h_values
}

# --- Spatial block bootstrap ---
# Resamples CV folds (not individual observations) with replacement to
# preserve spatial structure. Each iteration draws n_folds fold indices,
# pools their training data, and computes the requested statistic.
manual_spatial_bootstrap <- function(data, response, predictors, params,
                                     spatial_folds, stat_function, R = 1000,
                                     additional_args = NULL) {
  n_folds <- length(spatial_folds)
  results <- list()
  
  for (i in seq_len(R)) {
    if (i %% 100 == 0) message("  Bootstrap iteration ", i, "/", R)
    
    sampled_folds <- sample(seq_len(n_folds), size = n_folds, replace = TRUE)
    train_indices <- unlist(lapply(sampled_folds, function(idx) spatial_folds[[idx]][[1]]))
    bootstrap_data <- data[train_indices, ]
    
    if (nrow(bootstrap_data) < 10) next
    
    tryCatch({
      args <- c(list(data = bootstrap_data, response = response,
                     predictors = predictors, params = params),
                additional_args)
      results[[i]] <- do.call(stat_function, args)
    }, error = function(e) {
      message("  Iteration ", i, " failed: ", e$message)
    })
  }
  
  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) == 0) stop("All bootstrap iterations failed")
  
  if (is.list(results[[1]])) return(results)
  if (is.matrix(results[[1]]) || is.data.frame(results[[1]])) return(results)
  do.call(rbind, results)
}
