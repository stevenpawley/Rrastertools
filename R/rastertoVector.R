#' Raster to vector conversion with generalization
#'
#' @param rasterlayer RasterLayer object
#' @param gapp character, path to GRASS GIS
#' @param rmarea_km2 Minimum size of polygons in km2
#' @param simplify_thres Neuuman reduction simplification threshold
#' @param smoothing_iter Number of chaiken smoothing iterations
#' @param cycles Number of reduction and smoothing iterations
#'
#' @return
#' @export
rastertoVector = function(rasterlayer, gapp = NULL, rmarea_km2, simplify_thres, smoothing_iter, cycles) {
  
  # initiate a temporary grass session
  if (missing(gapp)) {
    link2GI::linkGRASS7(x = rasterlayer)
  } else {
    rgrass7::initGRASS(gisBase = gapp, override = T)  
  }
  
  # convert raster to vector and remove small areas
  outfile = tempfile(fileext = '.tif')
  raster::writeRaster(x = rasterlayer, filename = outfile)
  rgrass7::execGRASS('r.in.gdal', input=outfile, output='rasterlayer', flags='o')
  rgrass7::execGRASS('g.region', raster='rasterlayer')
  rgrass7::execGRASS('r.to.vect', flags=c('s', 'overwrite'), input='rasterlayer',
                     type='area', output='rasterlayer')
  rgrass7::execGRASS('v.clean', input='rasterlayer', output='rasterlayerclean',
                     tool='rmarea', threshold=rmarea_km2*1000000, flag='overwrite')
  
  # repeatedly simplify and smooth
  for (i in 1:cycles) {
    invector = 'rasterlayerclean'
    rgrass7::execGRASS('v.generalize', input=invector, output='rasterlayersimple',
                       method='reduction', threshold=simplify_thres, flags='overwrite')
    rgrass7::execGRASS('v.generalize', input='rasterlayersimple', output=invector,
                       method='chaiken', threshold=1000000,
                       iterations=smoothing_iter, flags='overwrite')
  }
  
  # remove small areas
  rgrass7::execGRASS('v.clean', input=invector, output='final',
                     tool='rmarea', threshold=rmarea_km2*1000000, flag='overwrite')
  
  generalized_py = rgrass7::readVECT(vname = 'final') %>% sf::st_as_sf()
  
  return(generalized_py)
}
