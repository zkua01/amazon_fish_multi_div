#extract-stats
calc_stats <- function(monthly_df) {
  #  monthly_values <- exact_extract(raster_stack, basins, fun = "mean")
  #  monthly_values_df <- as.data.frame(monthly_values)
  annual_means <- rowMeans(monthly_df, na.rm = TRUE)
  annual_sums <- rowSums(monthly_df, na.rm = TRUE)
  max_values <- apply(monthly_df, 1, max, na.rm = TRUE)
  min_values <- apply(monthly_df, 1, min, na.rm = TRUE)
  cv_values <- apply(monthly_df, 1, function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  })
  
  data.frame(
    mean = annual_means,
    annual = annual_sums,
    max = max_values,
    min = min_values,
    cv = cv_values
  )
}
