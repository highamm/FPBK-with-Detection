#' Summarizes a fitted spatial linear model
#'
#' @param object is an object generated from \code{\link{slmfit}}.
#' @param ... are options to be passed to \code{print}
#' @return a list with \itemize{
#'   \item model formula
#'   \item a table of fixed effects estimates and associated standard errors. For models with detection, these are the fixed effects estimates for the abundance model. Due to the complexity in estimating the fixed effects, standard errors and p-values are not given in the summary output when using detection from radiocollar data.
#'   \item estimated spatial covariance parameter estimates
#'   \item residuals
#'        }
#' @import stats
#' @export

summary.slmfit <- function(object, ...) {
  
  catcall <- object$FPBKpredobj$formula
  
  predictornames <- object$PredictorNames
  NAvec <- rep(NA, times = length(predictornames))
  
  regcoefs <- object$CoefficientEsts
  regvar <- object$BetaCov
  p <- length(regcoefs)
  n <- object$SampSize
  
  sereg <- sqrt(diag(as.matrix(regvar)))
  
  tvec <- NAvec
  tvec <- regcoefs / sereg
  pvec <- NAvec
  pvec <- round(100000 * (1 - pt(abs(regcoefs / sereg),
    df = n - p)) * 2) / 100000
  
  fixed.eff.est <- data.frame(##FactorLevel = predictornames,
    Estimate = regcoefs,
    std.err = sereg, t.value = tvec, prob.t = pvec)
  fixed.effects.estimates = fixed.eff.est
  
  covmodels <- object$SpatialParmEsts
  covmodelout <- data.frame(covmodels)
  colnames(covmodelout) <- paste(object$CovarianceMod, "Model")
  
  resid_vec <- object$resids
  ##residualsum <- c(min(residuals), quantile(residuals, c(0.25, 0.5,
  ##  0.75)), max(residuals))
  ##generalizedr2 <- GR2(object)
  
  outpt <- list(catcall = catcall,
    fixed.effects.estimates = fixed.effects.estimates,
    covariance.parameters = covmodelout,
    resid_vec)
    ##generalizedr2)
  names(outpt) <- c("catCall", "FixedEffects", "CovarianceParms",
    "Residuals") ##"GeneralizedR2")
  class(outpt) <- "summary.slmfit"
  return(outpt)
  
}
