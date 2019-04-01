#' Performs Block Kriging for Multiple Strata Separately
#'
#' Runs \code{slmfit}, \code{predict}, and \code{FPBKoutput} for each
#' stratum in a user-specified stratification variable. Note: if fitting the spatial model separately for each stratum, then stratum cannot be included as a covariate in the \code{formula} argument. 
#'
#' @param formula is an R linear model formula specifying density as the
#' response variable as well as covariates for predicting densities on the unsampled sites.
#' @param data is the data set with the response column of densities, the covariates to
#' be used for the block kriging, and the spatial coordinates for all of the sites.
#' @param xcoordcol is the name of the column in the data frame with x coordinates or longitudinal coordinates
#' @param ycoordcol is the name of the column in the data frame with y coordinates or latitudinal coordinates
#' @param CorModel is the covariance structure. By default, \code{CorModel} is
#' Exponential but other options include the Spherical and Gaussian.
#' @param coordtype specifies whether spatial coordinates are in latitude, longitude (\code{LatLon}) form or UTM (\code{UTM}) form.
#' @param estmethod is either the default \code{"REML"} for restricted
#' maximum likelihood to estimate the covariance parameters and
#' regression coefficients or \code{"ML"} to estimate the covariance
#' parameters and regression coefficients.
#' @param covestimates is an optional vector of covariance parameter estimates (nugget, partial sill, range). If these are given and \code{estmethod = "None"}, the the provided vector are treated as the estimators to create the covariance structure.
#' @param detectionobj is a fitted model obj from \code{get_detection}. The default is for this object to be \code{NULL}, resulting in
#' spatial prediction that assumes perfect detection.
#' @param areacol is the name of the column with the areas of the sites. By default, we assume that all sites have equal area, in which
#' case a vector of 1's is used as the areas.
#' @param stratcol is the column in the data set that contains the stratification variable. 
#' @return a report with information about the predicted total across all sites as well as variogram information for each stratum in the column \code{stratcol}.
#' @import stats
#' @export multistrat


multistrat <- function(formula, data, xcoordcol, ycoordcol,
  CorModel = "Exponential",
  coordtype = "LatLon", estmethod = "REML",
  covestimates = c(NA, NA, NA),
  detectionobj = NULL,
  areacol = NULL,
  stratcol = NULL) {
  
  data$stratvar <- factor(data[ ,stratcol])

  slmfitouts <- vector("list",  length(levels(data$stratvar)))
  predictouts <- vector("list",  length(levels(data$stratvar)))
  dfout <- vector("list", length(levels(data$stratvar)))
  prediction <- vector("list", length(levels(data$stratvar)))
  predictionvar <- vector("list", length(levels(data$stratvar)))
  
  for (k in 1:length(levels(data$stratvar))) {
  slmfitouts[[k]] <- slmfit(formula = formula,
    data = subset(data, data$stratvar == levels(data$stratvar)[k]),
    xcoordcol = xcoordcol, ycoordcol = ycoordcol, CorModel = CorModel, 
    coordtype = coordtype, estmethod = estmethod,
    covestimates = covestimates, detectionobj = detectionobj,
    areacol = areacol)
  
  predictouts[[k]] <- predict(object = slmfitouts[[k]],
    FPBKcol = NULL, detinfo = c(1, 0))
  
  ## essentially have to store the report for each k,
  ## then append these reports to the final report for all
  ## of the strata that is generated outside of the loop
  FPBKoutput(predictouts[[k]], conf_level = c(0.80, 
    0.90, 0.95),
    get_krigmap = FALSE, get_sampdetails = TRUE,
    get_variogram = TRUE, get_report = TRUE,
    nbreaks = 4,
    breakMethod = 'quantile', 
    pointsize = 2)
  
  prediction[[k]] <- predictouts[[k]]$FPBK_Prediction
  predictionvar[[k]] <- predictouts[[k]]$PredVar
  dfout[[k]] <- predictouts[[k]]$Pred_df
  ## could have a page for the low stratum report, a page for the
  ## high stratum report, and then a combined report....it will take 
  ## a bit of work combining the results but not too much I think
  }
  
  allFPBKobj <- vector("list", length(predictouts[[k]]))
  allFPBKobj[[1]] <- matrix(do.call("sum", prediction))
  ## need to update the variance to account for covariance
  ## from using the same detection for each stratum
  allFPBKobj[[2]] <- matrix(do.call("sum", predictionvar))
  allFPBKobj[[3]] <- do.call("rbind", dfout)
  allFPBKobj[[4]] <- c(NA, NA, NA)
  allFPBKobj[[5]] <- formula
  allFPBKobj[[6]] <- "DontNeed"

  names(allFPBKobj) <- c("FPBK_Prediction", "PredVar",
    "Pred_df", "SpatialParms", "formula", "CorModel")
  
  FPBKoutput(allFPBKobj, conf_level = c(0.80, 
    0.90, 0.95),
    get_krigmap = TRUE, get_sampdetails = TRUE,
    get_variogram = FALSE, get_report = TRUE,
    nbreaks = 4,
    breakMethod = 'quantile', 
    pointsize = 7)

  
 ## ggplot2::scale_shape_manual("Samp Indicator", 
  ##    labels = c("Unsampled", "Sampled", "x", "xx"),
  ##    values = shapevals)

}


