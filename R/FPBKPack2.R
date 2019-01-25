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
#' \code{FPBKPack2} Main Functions:
#' 
#' \code{FPBKpred} is the function used to generate the kriged predictions
#' on the unsampled sites.
#' 
#' \code{FPBKoutput} takes the resulting object from FPBKpred to construct
#' confidence intervals, variograms, and maps.
#' 
#' @docType package
#' @name FPBKPack2
NULL