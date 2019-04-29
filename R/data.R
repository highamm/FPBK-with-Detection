#' Data Set with Uncorrelated Poisson Counts.
#'
#' A data set that can be used with the \code{FPBKPack2} package. In 
#' this example, the counts are uncorrelated, the covariates are
#' generated as uniform random variables, and the sites fall on a 
#' regular grid.
#'
#' @format A data frame with 40  rows and 7 variables:
#' \describe{
#'   \item{counts}{counts, with NA values for unsampled sites}
#'   \item{pred1}{a possible predictor}
#'   \item{pred2}{a second possible predictor}
#'   \item{xcoords}{coordinates on the x-axis}
#'   \item{ycoords}{coordinates on the y-axis}
#'   \item{dummyvar}{an extra variable}
#'   \item{areavar}{Variable for the area of each plot}
#'   ...
#' }
"exampledataset"


#' Data Set with Alaska Moose Observations.
#'
#' A data set that can be used with the \code{sptotal} package. In
#' this example, the counts are of moose on 860 sites of equal area.
#'
#' @format A spatial polygons object inclding:
#' \describe{
#'   \item{CENTRLAT}{The latitude of the centroid for each site}
#'   \item{CENTRLON}{The latitude of the centroid for each site}
#'   \item{STRAT}{A stratification variable}
#'   \item{TOTAL}{The total moose count on each site}
#'   ...
#' }
"AKmoose"