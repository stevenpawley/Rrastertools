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
#' @param sf_obj sf object
#' @param field name of attribute that defines categories in the geomap to
#' calculate grid distances to
#' @param rasterlayer RasterLayer object to use as a template for the grid
#' distances
#' @param geomap n_jobs, number of processing cores, default = 1
#'
#' @return RasterStack of grid distances
#' @export
proximityMapFeatures <- function(sf_obj, field, rasterlayer, n_jobs = 1) {

  # some checks
  if (is(rasterlayer, "RasterStack") | is(rasterlayer, "RasterBrick")) {
    rasterlayer <- rasterlayer[[1]]
  }

  if (is(sf_obj, "sf") == FALSE) {
    stop("geomap must be a simple features dataframe")
  }

  if (missing(field)) {
    stop("field must be specified")
  }

  # rasterize shapes based on field attribute
  rasterized_shapes <- fasterize::fasterize(
    raster = rasterlayer,
    sf = sf_obj,
    by = field
  )

  raster_stats <- sapply(1:raster::nlayers(rasterized_shapes), function(i)
    raster::maxValue(rasterized_shapes[[i]]))
  invalid_layers <- which(is.na(raster_stats))

  rasterized_shapes <- raster::dropLayer(rasterized_shapes, invalid_layers)

  # calculate distances
  dist_fun <- function(i) {
    prox <- raster::distance(rasterized_shapes[[i]])
    names(prox) <- paste("proximity", names(rasterized_shapes)[[i]], sep = ".")
    return(prox)
  }

  if (n_jobs == 1) {
    proximities <- lapply(1:raster::nlayers(rasterized_shapes), dist_fun)
  } else {
    cl <- parallel::makeCluster(n_jobs)
    clusterExport(cl, "rasterized_shapes")
    clusterEvalQ(cl, "raster")
    proximities <- parallel::parLapply(cl, 1:raster::nlayers(rasterized_shapes), dist_fun)
    parallel::stopCluster(cl)
  }
  proximities <- raster::stack(proximities)

  return(proximities)
}

#' Corner-based Euclidean Distance Fields
#' 
#' Calculates corner-based related euclidean distance fields
#' (Behrens et al., 2018) to use as predictors in spatial models
#'
#' @param object RasterLayer object to use as a template
#'
#' @return RasterStack containing corner and centre coordinate EDM grids
#' @export
#'
#' @examples
euclideanDistanceFields <- function(object) {
  
  ext <- raster::extent(object)
  
  # top left
  topleft <- raster::raster(nrows = raster::nrow(object),
                            ncols = raster::ncol(object),
                            resolution = raster::res(object),
                            crs = raster::crs(object),
                            ext = ext)
  topleft[1, 1] <- 1
  topleft <- raster::distance(topleft)
  
  # top right
  topright <- raster::raster(nrows = raster::nrow(object),
                             ncols = raster::ncol(object),
                             resolution = raster::res(object),
                             crs = raster::crs(object),
                             ext = ext)
  topright[1, ncol(object)] <- 1
  topright <- raster::distance(topright)
  
  # bottom left
  bottomleft <- raster::raster(nrows = raster::nrow(object),
                               ncols = raster::ncol(object),
                               resolution = raster::res(object),
                               crs = raster::crs(object),
                               ext = ext)
  bottomleft[nrow(object), 1] <- 1
  bottomleft <- raster::distance(bottomleft)
  
  # bottom right
  bottomright <- raster::raster(nrows = raster::nrow(object),
                                ncols = raster::ncol(object),
                                resolution = raster::res(object),
                                crs = raster::crs(object),
                                ext = ext)
  bottomright[nrow(object), ncol(object)] <- 1
  bottomright <- raster::distance(bottomright)
  
  # centre
  centre <- raster::raster(nrows = raster::nrow(object),
                           ncols = raster::ncol(object),
                           resolution = raster::res(object),
                           crs = raster::crs(object),
                           ext = ext)
  centre[as.integer(nrow(object)/2), as.integer(ncol(object)/2)] <- 1
  centre <- raster::distance(centre)
  
  EDM <- raster::stack(topleft, topright, centre, bottomleft, bottomright)
  EDM <- setNames(EDM,
                  c('EDM_topleft',
                    'EDM_topright',
                    'EDM_centre',
                    'EDM_bottomleft',
                    'EDM_bottomright'))
  return(EDM)
}


#' Produces RasterLayer objects filled with rotated coordinate values
#'
#' @param object RasterLayer object
#' @param n_angles numeric, number of angles to rotate coordinates by
#'
#' @return RasterStack object
#' @export
rotatedCoordinateGrids <- function(object, n_angles) {
  anglegrids <- as(object, "SpatialGridDataFrame")
  angles <- NISTunits::NISTdegTOradian(seq(from = 0, to = 180, length.out = n_angles))

  for (i in seq_along(angles)) {
    newlayer <- paste0("angle", i)
    anglegrids[[newlayer]] <- coordinates(object)[, 1] + angles[i] * coordinates(object)[, 2]
  }

  anglegrids <- anglegrids[2:ncol(anglegrids)]
  anglegrids <- raster::stack(anglegrids)

  return(anglegrids)
}


#' Creates RasterLayer objects filled by the x and y grid coordinates
#'
#' @param object RasterLayer object
#'
#' @return RasterStack object
#' @export
xyCoordinateGrids <- function(object) {
  xy_coords <- raster::xyFromCell(object, cell = 1:raster::ncell(object))
  object <- raster::stack(object)

  object[["xgrid"]] <- raster(
    nrows = nrow(object),
    ncols = ncol(object),
    ext = extent(object),
    crs = crs(object),
    vals = xy_coords[, 1]
  )

  object[["ygrid"]] <- raster(
    nrows = nrow(object),
    ncols = ncol(object),
    ext = extent(object),
    crs = crs(object),
    vals = xy_coords[, 2]
  )

  object <- object[[c("xgrid", "ygrid")]]

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
    ext <- raster::extent(y)
    xcells <- raster::ncol(y)
    ycells <- raster::nrow(y)

    # get dimensions to perform the estimation over from bbox of data.points
  } else {
    ext <- sp::bbox(data.points)

    if (is.null(xcells) | is.null(ycells)) {
      stop("Need to supply xcells and ycells if a raster object is not supplied")
    }
  }

  xmin <- ext[1, 1]
  xmax <- ext[1, 2]
  ymin <- ext[2, 1]
  ymax <- ext[2, 2]

  # compute the 2D binned kernel density estimate
  est <- KernSmooth::bkde2D(
    coords,
    bandwidth = c(selected.xbandwidth, selected.ybandwidth),
    gridsize = c(xcells, ycells),
    range.x = list(
      c(xmin, xmax),
      c(ymin, ymax)
    )
  )

  # create raster
  est.raster <- raster::raster(
    list(
      x = est$x1,
      y = est$x2,
      z = est$fhat
    )
  )
  raster::projection(est.raster) <- sp::crs(data.points)

  return(est.raster)
}


#' Creates a RasterStack object containing euclidean distances to points
#'
#' @param points sf object
#' @param rasterlayer RasterLayer object to use as grid template
#' @param field character, name of field to group points into classes
#' @param n_classes numeric, number of classes
#' @param method character, either 'equal_intervals' or 'quantiles'
#' Used to create groups of points based on the field attribute
#'
#' @return RasterStack object of point distances
#' @export
distanceFromPointIntervals <- function(points, rasterlayer, field = NULL,
                                       n_classes = 10,
                                       method = "equal_intervals") {
  if (missing(field)) {
    stop("Field attribute must be supplied")
  }

  # create point buffers
  if (method == "equal_intervals") {
    classes <- cut(points[[field]],
      breaks = seq(min(points[[field]]),
        max(points[[field]]),
        length = n_classes
      )
    )
  } else if (method == "quantiles") {
    classes <- cut(points[[field]],
      breaks = stats::quantile(points[[field]],
        probs = seq(0, 1, length.out = n_classes),
        include.lowest = TRUE
      )
    )
  }

  buffer_grids <- lapply(
    split(picks, classes),
    function(points, raster_grid) {
      raster::distanceFromPoints(raster_grid, xy = as(points, "Spatial"))
    },
    raster_grid = rasterlayer
  ) %>% stack()

  names(buffer_grids) <- paste0("buffer", seq(1, nlayers(buffer_grids)))

  return(buffer_grids)
}


#' Sample-based euclidean distance fields
#' 
#' Calculates sample-based euclidean distance fields, i.e. buffer distances
#' to each point in a simple features object
#'
#' @param object RasterLayer to use as template
#' @param sf_obj Simple features object containing POINT geometries
#' @param n_jobs numeric, optionally use parallel calculation over multiple cores
#'
#' @return RasterStack of sample-based EDMs
#' @export
#'
#' @examples
sampleEuclideanDistanceFields <- function(object, sf_obj, n_jobs = 1) {
  
  if (n_jobs == 1) {
    
    buffer_grids <- lapply(seq_along(sf_obj), function(i)
      raster::distanceFromPoints(object, xy = as(sf_obj[i, ], "Spatial")))
    
  } else {
    
    cl <- parallel::makeCluster(n_jobs)
    parallel::clusterEvalQ(cl, "raster")
    parallel::clusterExport(cl, c("object", "sf_obj"))

    sf_obj <- as(sf_obj, "Spatial")
    
    buffer_grids <- parallel::parLapply(cl, seq_along(sf_obj), function(i)
      raster::distanceFromPoints(object, xy = sf_obj[i, ]))
    
    parallel::stopCluster(cl)
  }
  
  buffer_grids <- stack(buffer_grids)
  buffer_grids <- setNames(buffer_grids, paste0("buffer", seq(nlayers(buffer_grids))))
  
  return(buffer_grids)
}

