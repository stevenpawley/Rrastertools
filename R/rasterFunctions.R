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
#' @param parallel logical, use parallel processing. Default is TRUE
#' @param n_jobs numeric, number of processor cores. Default is ncores-1
#'
#' @return RasterStack of grid distances
#' @export
proximityMapFeatures = function(geomap, field, rasterlayer, parallel = TRUE, 
                                n_jobs = NULL) {
  
  f = function(X, geomap, field, rasterlayer, ...) {
    # function to calculate proximity to simple feature geometries
    
    selected_poly = geomap[geomap[[field]] == X, ]
    
    if (nrow(selected_poly) > 0) {
      ras = raster::rasterize(selected_poly, rasterlayer)
      dist_to_ras = raster::distance(ras, doEdge = TRUE)
    } else {
      dist_to_ras = NA
    }
    
    return(dist_to_ras)
  }
  
  # get the categories to split the feature by
  classes = as.character(unique(geomap[[field]]))
  
  if (parallel == TRUE) {
  
    tryCatch(
      expr = {
        # initiate a cluster with the required packages
        if (missing(n_jobs))
          n_jobs = parallel::detectCores()-1
        
        cl = parallel::makeCluster(n_jobs)
        clusterEvalQ(cl, library(sf))
        clusterEvalQ(cl, library(raster))

        # calculate proximity to features in parallel
        proximities = parallel::parLapply(
          cl = cl,
          X = classes,
          fun = f, geomap=geomap, field=field, rasterlayer=rasterlayer)
        },
      
      finally = {
        # close the cluster
        parallel::stopCluster(cl)
        }
      )

    
    # single-threaded version
  } else {
    proximities = lapply(
      X = classes,
      FUN = f, geomap, field, rasterlayer)
  }
  
  proximities = stack(proxmities)
  proximities = setNames(proximities, classes)
  
  return(stack(proximities))
}


#' Produces RasterLayer objects filled with rotated coordinate values
#'
#' @param object RasterLayer object
#' @param n_angles numeric, number of angles to rotate coordinates by
#'
#' @return RasterStack object
#' @export
rotatedCoordinateGrids = function(object, n_angles) {
  
  anglegrids = as(object, "SpatialGridDataFrame")
  angles = NISTunits::NISTdegTOradian(seq(from = 0, to = 180, length.out = n_angles))
  
  for (i in seq_along(angles)) {
    newlayer = paste0('angle', i)
    anglegrids[[newlayer]] = coordinates(object)[,1] + angles[i] * coordinates(object)[,2]
  }
  
  anglegrids = anglegrids[2:ncol(anglegrids)]
}


#' Creates RasterLayer objects filled by the x and y grid coordinates
#'
#' @param object RasterLayer object
#'
#' @return RasterStack object
#' @export
xyCoordinateGrids = function(object) {
  
  xy_coords = raster::xyFromCell(object, cell = 1:raster::ncell(object))
  
  object[['xgrid']] = raster(
    nrows = nrow(gr),
    ncols = ncol(gr),
    ext = extent(gr),
    crs = crs(gr),
    vals = xy_coords[,1])
  
  object[['ygrid']] = raster(
    nrows = nrow(gr),
    ncols = ncol(gr),
    ext = extent(gr),
    crs = crs(gr),
    vals = xy_coords[,2])
  
  return(object)
}


#' 2d kernel density estimation for raster data
#'
#' @param data.points SpatialPointsDataFrame object
#' @param y RasterLayer object to use as a template for the output of the
#' kernel density estimator, optional.
#' @param xcells number of grid cells in x dimension of output raster
#' @param ycells number of grid cells in y dimension of output raster
#' @return
#' @export
kernelDensity2D <- function(data.points, y = NULL, xcells = NULL, ycells = NULL) {
  
  # get the coordinates
  coords <- sp::coordinates(data.points)
  
  # bandwidth selection
  selected.xbandwidth <- KernSmooth::dpik(coords[, 1])
  selected.ybandwidth <- KernSmooth::dpik(coords[, 2])
  
  # get dimensions to perform the estimation over from raster if supplied
  if (!is.null(y)) {
    ext = raster::extent(y)
    xcells = raster::ncol(y)
    ycells = raster::nrow(y)
    
    # get dimensions to perform the estimation over from bbox of data.points
  } else {
    ext = sp::bbox(data.points)
    
    if (is.null(xcells) | is.null(ycells))
      stop('Need to supply xcells and ycells if a raster object is not supplied')
  }
  
  xmin = ext[1, 1]
  xmax = ext[1, 2]
  ymin = ext[2, 1]
  ymax = ext[2, 2]
  
  # compute the 2D binned kernel density estimate
  est <- KernSmooth::bkde2D(
    coords,
    bandwidth = c(selected.xbandwidth, selected.ybandwidth),
    gridsize = c(xcells, ycells),
    range.x = list(c(xmin, xmax),
                   c(ymin, ymax))
  )
  
  # create raster
  est.raster = raster::raster(
    list(
      x = est$x1,
      y = est$x2,
      z = est$fhat
    ))
  raster::projection(est.raster) <- sp::crs(data.points)
  
  return (est.raster)
}


#' Hillshade PCA function
#'
#' @param x RasterLayer object representing a DEM
#' @param azi numeric vector of azimuth of different shading directions
#' @param n numeric, number of pixels to sample for PCA. Default = 5000
#'
#' @return list, containing hillshade pca and the associated prcomp class
#' @export

hillshadePCA = function(x, azi = seq(0, 180, 22.5)[1:8], n = 5000) {
  # generate hillshades
  slope = raster::terrain(x, opt = 'slope')
  aspect = raster::terrain(x, opt = 'aspect')
  azi_shades = stack(lapply(azi, function(x) raster::hillShade(
    slope, aspect, angle=x, normalize = T)))
  
  # train pca
  sr = raster::sampleRandom(azi_shades, n)
  pca = stats::prcomp(sr) 
  
  # create new rasters from PCA predictions
  hillshade_pca = stats::predict(azi_shades, pca, index=1:3)
  
  return(list(hillshade=hillshade_pca, pca=pca))
}
