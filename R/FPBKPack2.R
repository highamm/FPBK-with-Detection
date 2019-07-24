#' FPBKPack2: A package used for performing finite population block
#' kriging on polygonal count data. 
#' 
#' The package provides an option to perform FPBK on counts assuming
#' perfect detection or to perform an adjusted FPBK procedure to
#' correct for the possibility of imperfect detection. Imperfect
#' detection estimates can arise from a survey of radiocollared
#' survey units or, alternatively, from another method that gives
#' the estimated mean detection probability with a standard error
#' of that estimator.
#' 
#' \code{FPBKPack2} Main Functions
#' \describe{
#' \item{\code{\link{slmfit}}}{fits a spatial linear model to the observed data,}
#' \item{\code{\link{predict}}}{uses Finite Population Block Kriging to predict counts on unobserved sites,}
#' \item{\code{\link{FPBKoutput}}}{takes the resulting object from FPBKpred to construct confidence intervals, variograms, and maps,}
#' \item{\code{\link{get_reportdoc}}}{puts the output into nice HTML format, and}
#' \item{\code{\link{multistrat}}}{combines \code{\link{slmfit}}, \code{\link{predict}}, and \code{\link{FPBKoutput}} into one function so that the user can fit the spatial linear model to strata separately.}
#' }
#' @docType package
#' @name FPBKPack2
NULL