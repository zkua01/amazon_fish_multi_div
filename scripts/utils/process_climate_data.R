#modified: process rasters sequentially and extract stats incrementally
process_climate_data <- function(variable_name, file_path, basins_object, file_structure = "flat"){
  #initialize empty datadframe for results
  monthly_results <- data.frame(matrix(nrow = nrow(basins_object), ncol = 0))
  basins_extent <- terra::ext(basins_raw)+0.1
  
  #process each monthly file sequentially
  for(month in months){
    message("Processing ", variable_name, " month ", month, "...")
    
    if (file_structure == "nested"){
      folder_path <- file.path(file_path, paste0(variable_name,"_",as.numeric(month)))
      file <- file.path(folder_path, "w001001.adf")
    } else {
      file <- paste0(file_path, "/", variable_name, "_", month, ".tif")
    }
    message("Reading: ", file)
    
    if (!file.exists(file)){
      stop("File not found: ", file)
    }
    #read and process single monthly file
    
    r <- crop(rast(file), basins_extent)
    if (!is.na(crs(r)) && !grepl("EPSG:5880", crs(r))){
      message("Reprojecting cropped raster: ", basename(file))
      r <- project(r, "EPSG:5880", method = "bilinear")
    }
    
    message("Extracting values...")
    #extract mean value for this month
    month_values <- exact_extract(r, basins_object, fun = "mean")
    monthly_results[[paste0("month_",month)]] <- month_values
    
    rm(r, month_values)
    gc()
    
    message("Completed month ", month)
  }
  calc_stats(monthly_results)
}
