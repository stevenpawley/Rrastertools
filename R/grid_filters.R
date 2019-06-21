#' Morphological filter
#'
#' @param rasterlayer RasterLayer object
#' @param size Size of structuring element (in grid cells)
#' @param type Type of filter c('erosion', 'dilation', 'opening', 'closing)
#'        Currently only closing is supported
#' @param parallel Logical, run focal function in parallel (default = FALSE)
#'
#' @return Rasterlayer after morphological filtering
#' @export
#' @importFrom raster focal
#' @importFrom spatial.tools sfQuickInit rasterEngine sfQuickStop
#' @importFrom parallel detectCores
morphological_filter <- function(rasterlayer, size = 3, type = "closing", parallel = F) {

  # single core morphological closing using the raster package
  if (parallel == F) {
    f <- matrix(1, size, size)

    if (type == "closing") {
      dilated <- focal(rasterlayer, w = f, fun = max)
      filtered <- focal(dilated, w = f, fun = min)
    }


    # parallelized morphological closing using spatial.tools package
  } else {
    if (type == "closing") {
      dilation <- function(inraster) maxfocal <- apply(inraster, 3, max)
      erosion <- function(inraster) maxfocal <- apply(inraster, 3, min)

      sfQuickInit(cpus = detectCores())
      dilated <- rasterEngine(
        inraster = rasterlayer, fun = dilation,
        window_dims = c(size, size)
      )
      filtered <- rasterEngine(
        inraster = dilated, fun = erosion,
        window_dims = c(size, size)
      )
      sfQuickStop()
    }
  }

  return(filtered)
}
