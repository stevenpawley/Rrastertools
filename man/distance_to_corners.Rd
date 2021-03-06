% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grid_distances.R
\name{distance_to_corners}
\alias{distance_to_corners}
\title{Corner-based Euclidean Distance Fields}
\usage{
distance_to_corners(object)
}
\arguments{
\item{object}{RasterLayer object to use as a template}
}
\value{
RasterStack containing corner and centre coordinate EDM grids
}
\description{
Calculates corner-based related euclidean distance fields
(Behrens et al., 2018) to use as predictors in spatial models
}
