#' Corner-based Euclidean Distance Fields
#' 
#' Calculates corner-based related euclidean distance fields
#' (Behrens et al., 2018) to use as predictors in spatial models
#'
#' @param object RasterLayer object to use as a template
#'
#' @return RasterStack containing corner and centre coordinate EDM grids
#' @export
#' @importFrom raster raster res crs distance stack
#' @importFrom stats setNames
distance_to_corners <- function(object) {
  
  ext <- raster::extent(object)
  
  # top left
  topleft <- raster(
    nrows = nrow(object),
    ncols = ncol(object),
    resolution = res(object),
    crs = crs(object),
    ext = ext)
  topleft[1, 1] <- 1
  topleft <- distance(topleft)
  
  # top right
  topright <- raster(
    nrows = nrow(object),
    ncols = ncol(object),
    resolution = res(object),
    crs = crs(object),
    ext = ext)
  topright[1, ncol(object)] <- 1
  topright <- distance(topright)
  
  # bottom left
  bottomleft <- raster(
    nrows = nrow(object),
    ncols = ncol(object),
    resolution = res(object),
    crs = crs(object),
    ext = ext)
  bottomleft[nrow(object), 1] <- 1
  bottomleft <- distance(bottomleft)
  
  # bottom right
  bottomright <- raster(
    nrows = nrow(object),
    ncols = ncol(object),
    resolution = res(object),
    crs = crs(object),
    ext = ext)
  bottomright[nrow(object), ncol(object)] <- 1
  bottomright <- distance(bottomright)
  
  # centre
  centre <- raster(
    nrows = nrow(object),
    ncols = ncol(object),
    resolution = res(object),
    crs = crs(object),
    ext = ext)
  centre[as.integer(nrow(object)/2), as.integer(ncol(object)/2)] <- 1
  centre <- distance(centre)
  
  EDM <- stack(topleft, topright, centre, bottomleft, bottomright)
  EDM <- setNames(EDM,
    c('edm_topleft', 'edm_topright', 'edm_centre', 'edm_bottomleft', 'edm_bottomright'))
  return(EDM)
}


#' Creates RasterLayer objects filled by the x and y grid coordinates
#'
#' @param object RasterLayer object
#'
#' @return RasterStack object
#' @export
#' @importFrom methods is
#' @importFrom raster stack xyFromCell ncell raster extent crs
coordinate_grids <- function(object) {
  
  if (is(object, "RasterStack") | is(object, "RasterBrick"))
    object <- object[[1]]
  
  xy_coords <- xyFromCell(object, cell = 1:ncell(object))
  object <- stack(object)

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


#' Sample-based euclidean distance fields
#' 
#' Calculates sample-based euclidean distance fields, i.e. buffer distances
#' to each point in a simple features object
#'
#' @param object RasterLayer to use as template
#' @param sf_obj Simple features object containing POINT geometries
#'
#' @return RasterStack of sample-based EDMs
#' @export
#' @importFrom raster distanceFromPoints stack nlayers
#' @importFrom methods as
#' @importFrom stats setNames
distance_to_features <- function(object, sf_obj) {
  
  buffer_grids <- lapply(1:nrow(sf_obj), function(i)
    distanceFromPoints(object, xy = as(sf_obj[i, ], "Spatial")))

  buffer_grids <- stack(buffer_grids)
  buffer_grids <- setNames(buffer_grids, paste0("buffer", seq(nlayers(buffer_grids))))
  
  return(buffer_grids)
}
