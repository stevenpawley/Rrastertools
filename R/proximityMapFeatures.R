
#' Calculate grid distances to spatial features
#' 
#' This function splits a simple features object based on an attribute field
#' and rasterizes each category. Grid distances are then calculated to each 
#' categorical raster separately. This function represents an alternative to
#' one-hot encoding or using rasters that represent categorical features. 
#' However, instead of those features have distinct boundaries, the boundaries
#' here represent 'soft' margins defined by the proximity to the feature. This
#' avoids imprints of the polygons occurring in spatial predictions.
#'
#' @param geomap sf object
#' @param field name of attribute that defines categories in the geomap to 
#' calculate grid distances to
#' @param rasterlayer RasterLayer object to use as a template for the grid
#' distances
#'
#' @return RasterStack of grid distances
#' @export
proximityMapFeatures = function(geomap, field, rasterlayer, parallel = TRUE) {
  
  f = function(X, geomap, field, rasterlayer, ...) {
    selected_poly = geomap[geomap[[field]] == X, ]
    
    if (nrow(selected_poly) > 0) {
      ras = raster::rasterize(selected_poly, rasterlayer)
      dist_to_ras = raster::distance(ras, doEdge = TRUE)
    } else {
      dist_to_ras = NA
    }
    
    return(dist_to_ras)
  }
  
  classes = as.character(unique(geomap[[field]]))
  
  if (parallel == TRUE) {
  
    tryCatch(
      expr = {
        cl = parallel::makeCluster(getOption("cl.cores", parallel::detectCores()))
      
        proximities = parallel::parLapply(
          cl = cl,
          X = classes,
          fun = f, geomap=geomap, field=field, rasterlayer=rasterlayer)
        },
      
      finally = {
      parallel::stopCluster(cl)
        }
      )

    
  } else {
    proximities = lapply(
      X = classes,
      FUN = f, geomap, field, rasterlayer)
    
  }
  
  return(stack(proximities))
}
