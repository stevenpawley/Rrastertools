#' Merge factors containing a small number of observations
#'
#' Merges factors containing a small number of observations
#' in training data extracted from a raster
#' @param x RasterLayer object containing factors
#' @param y vector of extracted values to be used as training data
#' @param n integer, minimum number of observations required to retain a factor
#'
#' @return list
#' levels : data.frame with 3 columns containing raster values (ID),
#' the original levels that are associated with each value,
#' and  the collapsed levels
#' @export
merge_raster_factors <- function(x, y, n = 50) {
  
  # some checks
  if (is.null(levels(x)[[1]][, 2]))
    stop('RasterLayer does not contain any factor levels')
  
  if (ncol(levels(x)[[1]]) < 2 | names(levels(x)[[1]])[1] != 'ID')
    stop(paste('Raster attribute table needs to have at least two columns',
               'with ID providing the raster numeric pixel value',
               'and a second column providing the factor level'))
  
  if (is.null(levels(y)))
    stop('y does not contain any factor levels')
  
  # get list of all levels that exist in the RasterLayer
  lev <- levels(x)[[1]]
  
  # count the number of observations per level in y
  fact <- forcats::fct_count(y)
  lev["n"] <- fact["n"]
  
  if (n <= min(lev$n))
    stop(paste('None of the raster levels contain > n observations.',
               'This would cause all levels to be dropped. Increase n.'))
  
  # get factor levels to drop based on where n_obs < n
  levels_to_drop <- fact[which(fact$n < n), ][["f"]]
  
  # replace dropped raster levels with 'other'
  lev['collapsed'] <- lev[, 2]
  lev[, 'collapsed'] <- forcats::fct_other(
    f = lev[, 2], drop = levels_to_drop, other_level = "other")
  
  # remove the 'n' attributes
  lev <- lev[, which(names(lev) != 'n')]
  
  return(lev)
}
