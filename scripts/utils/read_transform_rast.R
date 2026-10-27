read_transform_rast <- function(path, target_crs = "EPSG:5880", method = "bilinear"){
  r <- rast(path)
  if(!is.na(crs(r)) && !grepl("5880", crs(r))) {
    message("Reprojecting raster: ", basename(path))
    r <- project(r, target_crs, method = method)
  }
  return(r)
}
