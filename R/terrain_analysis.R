#' Hillshade PCA function
#'
#' @param x RasterLayer object representing a DEM
#' @param azi numeric vector of azimuth of different shading directions
#' @param n numeric, number of pixels to sample for PCA. Default = 5000
#'
#' @return list, containing hillshade pca and the associated prcomp class
#' @export
#' @importFrom raster terrain stack hillShade sampleRandom
#' @importFrom stats prcomp predict
hillshade_pca <- function(x, azi = seq(0, 180, 22.5)[1:8], n = 5000) {
  # generate hillshades
  slope <- terrain(x, opt = "slope")
  aspect <- terrain(x, opt = "aspect")
  azi_shades <- stack(lapply(azi, function(x) hillShade(
    slope, aspect,
    angle = x, normalize = T
  )))
  
  # train pca
  sr <- sampleRandom(azi_shades, n)
  pca <- prcomp(sr)
  
  # create new rasters from PCA predictions
  hillshade_pca <- predict(azi_shades, pca, index = 1:3)
  
  return(list(hillshade = hillshade_pca, pca = pca))
}
