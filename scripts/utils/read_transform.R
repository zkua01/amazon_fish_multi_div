read_transform <- function(path, target_crs = 5880){
  obj <- st_read(path, quiet = T)
  if (st_crs(obj) != target_crs) {
    obj <- st_transform(obj, target_crs)
  }
  return(obj)
}
