#' Fits a Spatial Linear Model
#'
#' Estimates regression coefficients and spatial autocorrelation
#' parameters.
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
#' @return a list with \itemize{
#'   \item the spatial covariance estimates
#'   \item the regression coefficient estimates
#'   \item a list containing \enumerate{
#'        \item formula
#'        \item data
#'        \item xcoordcol
#'        \item ycoordcol
#'        \item CorModel
#'        \item Inverted covariance matrix on the sampled sites
#'        \item Covariance matrix on all sites
#'        }
#' }
#' @import stats
#' @export slmfit

slmfit <- function(formula, data, xcoordcol, ycoordcol,
  CorModel = "Exponential",
  coordtype = "LatLon", estmethod = "REML",
  covestimates = c(NA, NA, NA),
  detectionobj = NULL) {
  
  
  ## display error message if estmethod is set to None and the user
  ## does not input covariance estimates
  
  if (estmethod == "None" & sum(is.na(covestimates) > 0) > 0) {
    stop("If estmethod is set to None, then covestimates must
      be a vector without missing values
      with the estimated (nugget, partial sill, range)")
  }
  
  if (estmethod == "None" & length(covestimates) != 3) {
    stop("If estmethod is set to None, then covestimates must
      be a vector of length 3 with the estimated
      (nugget, partial sill, range)")
  }
  
  ## display some warnings if the user, for example tries to input the
  ## vector of xcoordinates as the input instead of the name of the column
  if(is.character(xcoordcol) == FALSE) {
    stop("xcoords must be a string giving the
      name of the column in the data set
      with the x coordinates")
  }
  if(is.character(ycoordcol) == FALSE) {
    stop("ycoords must be a string giving the
      name of the column in the data set
      with the y coordinates")
  }
  
  
  Xall <- model.matrix(formula, model.frame(formula, data,
    na.action = stats::na.pass))
  
  missingind <- base::apply(is.na(Xall), MARGIN = 1, FUN = sum)
  nmissing <- sum(missingind >= 1)
  
  datanomiss <- data[missingind == 0, ]
  
  ## give a warning if some of the predictors have missing values.
  if (nmissing >= 1) {
    warning(paste("There were", nmissing, "sites with predictors with missing values. These will be removed from the data set and further analysis will be completed without these observations."))
  }
  
  ## ASSUME that coordinates are lat/lon. Convert these to UTM
  if (coordtype != "LatLon" & coordtype != "UTM") {
    stop("coordtype must be a character string LatLon or UTM")
  } else if (coordtype == "LatLon") {
    xcoordsUTM <- LLtoUTM(cm = base::mean(datanomiss[ ,xcoordcol]),
      lat = datanomiss[ ,ycoordcol], lon = datanomiss[ ,xcoordcol])$xy[ ,1]
    ycoordsUTM <- LLtoUTM(cm = base::mean(datanomiss[ ,xcoordcol]),
      lat = datanomiss[ ,ycoordcol], lon = datanomiss[ ,xcoordcol])$xy[ ,2]
  } else if (coordtype == "UTM") {
    xcoordsUTM <- datanomiss[ ,xcoordcol]
    ycoordsUTM <- datanomiss[ ,ycoordcol]
  }
  
  detind <- is.null(detectionobj)
  
  if (is.null(detectionobj) == TRUE) {
  ## create the design matrix for unsampled sites, for all of the sites, and for the sampled sites, respectively.
  
  fullmf <- stats::model.frame(formula, na.action =
      stats::na.pass, data = datanomiss)
  yvar <- stats::model.response(fullmf, "numeric")
  density <- yvar
  
  ## remove any rows with missing values in any of the predictors
  formula.onlypreds <- formula[-2]
  
  X <- model.matrix(formula.onlypreds,
    model.frame(formula.onlypreds, datanomiss,
      na.action = stats::na.omit))
  
  ## divide data set into sampled sites and unsampled sites based
  ## on whether the response variable has a number (for sampled sites)
  ## or NA (for unsampled sites)
  
  ind.sa <- !is.na(yvar)
  ind.un <- is.na(yvar)
  data.sa <- datanomiss[ind.sa, ]
  data.un <- datanomiss[ind.un, ]
  
  m.un <- stats::model.frame(formula, data.un, na.action =
      stats::na.pass)
  
  
  Xu <- model.matrix(formula.onlypreds,
    model.frame(formula.onlypreds, data.un,
      na.action = stats::na.omit))
  
  
  ## sampled response values and design matrix
  m.sa <- stats::model.frame(formula, data.sa, na.action =
      stats::na.omit)
  z.sa <- stats::model.response(m.sa)
  Xs <- stats::model.matrix(formula, m.sa)
  z.density <- z.sa
  n <- nrow(Xs)
  
  
  prednames <- colnames(Xs)
  
  ## x and y coordinates for sampled and unsampled sites
  x.sa <- xcoordsUTM[ind.sa]
  y.sa <- ycoordsUTM[ind.sa]
  x.un <- xcoordsUTM[ind.un]
  y.un <- ycoordsUTM[ind.un]
  
  ## number of sites that were sampled
  n.sa <- nrow(Xs)
  ## number of sites that were not sampled
  n.un <- nrow(Xu)
  
  
  ## estimate the spatial parameters, the covariance matrix, and
  ## the inverse of the covariance matrix
  
  spat.est <- estcovparm(response = density,
    designmatrix = as.matrix(X),
    xcoordsvec = xcoordsUTM,
    ycoordsvec = ycoordsUTM, CorModel = CorModel,
    estmethod = estmethod,
    covestimates = covestimates)
  
  parms.est <- spat.est$parms.est
  Sigma <- spat.est$Sigma
  min2loglik <- spat.est$min2loglik
  
  nugget.effect <- parms.est[1]; parsil.effect <- parms.est[2]
  range.effect <- parms.est[3]
  
  
  Sigma.ssi <- solve(spat.est$qrV) / (nugget.effect + parsil.effect)
  
  ## the generalized least squares regression coefficient estimates
  
  betahat <- spat.est$b.hat
  
  ## estimator for the mean vector
  muhats <- Xs %*% betahat
  muhatu <- Xu %*% betahat
  
  resids <- z.sa - muhats
  
  muhat <- rep(NA, nrow(datanomiss))
  muhat[ind.sa == TRUE] <- muhats
  muhat[ind.sa == FALSE] <- muhatu
  
  
  ## returns a list with the following components:
  ## 1.) A vector of the estimated regression coefficients
  ## 2.) A vector of the estimated covariance parameters
  ## 3.) A list with the information needed by FPBKpred
  
  covparms <- as.vector(c(nugget.effect, parsil.effect, range.effect))
  betahatest <- as.vector(betahat)
  covest <- spat.est$covb
  
  
  names(covparms) <- c("Nugget", "Partial Sill", "Range")
  
  FPBKpredobj <- list(formula, datanomiss, xcoordsUTM, ycoordsUTM,
    estmethod, CorModel, Sigma, Sigma.ssi)
  names(FPBKpredobj) <- c("formula", "data", "xcoordsUTM",
    "ycoordsUTM",
    "estmethod","correlationmod", "covmat", "covmatsampi")
  obj <- list(covparms, betahatest, covest, min2loglik, prednames,
    n, CorModel, resids, Xs, z.sa, detind, FPBKpredobj)
  
  names(obj) <- c("SpatialParmEsts", "CoefficientEsts",
    "BetaCov", "minus2loglike", "PredictorNames", "SampSize",
    "CovarianceMod",
    "resids", "DesignMat", "Density", "DetectionInd",
    "FPBKpredobj")
  
  } else {
    ## fit with detection
    
    ## for now assume that the column names with the detection
    ## covariates are the same as those in get_detection
    
    ## estmethod must be maximum likelihood
    estmethod <- "ML"
    
    if (estmethod != "ML") {
      warning("If using detection, the estimation method must
        be Maximum Likelihood")
    }
    
    detectionform <- detectionobj$formula
    coefs <- detectionobj$coefs
    varmethod <- detectionobj$varmethod
    
    Xndetall <- model.matrix(detectionform[-2],
      model.frame(detectionform[-2], datanomiss,
      na.action = stats::na.pass))
    
    missingind <- base::apply(is.na(Xndetall), MARGIN = 1, FUN = sum)
    nmissing <- sum(missingind >= 1)
    
    ## data set containing complete observations
    ## for both the count predictors and the detection
    ## predictors
    datanomissboth <- datanomiss[missingind == 0, ]
    
    
    Xndet <- model.matrix(detectionform[-2],
      model.frame(detectionform[-2], datanomissboth,
        na.action = stats::na.pass))
    
    logitpi <- Xndet %*% coefs
    piest <- exp(logitpi) / (1 + exp(logitpi))
    
    if (detectionobj$varmethod == "Delta") {
      F <- detectionobj$covmat
    covmu <- Xndet %*% F %*% t(Xndet)
    
    vtv <- matrix(exp(logitpi) / (1 +
        exp(logitpi)) ^ 2, ncol = 1) %*%
      matrix(exp(logitpi) / (1 +
          exp(logitpi)) ^ 2, nrow = 1)
    
    Vnn <- vtv * covmu
    
    } else if (detectionobj$varmethod == "Bootstrap") {
      
    boot.coefs <- detectionobj$boot.coefs
    logitboot <- Xndet %*% t(boot.coefs)
    piboot <- exp(logitboot) / (1 + exp(logitboot))
    Vnn <- cov(t(piboot), use = "complete.obs")
    
    }
    
    
    fullmf <- stats::model.frame(formula, na.action =
        stats::na.pass, data = datanomissboth)
    yvar <- stats::model.response(fullmf, "numeric")
    density <- yvar
    
    ## remove any rows with missing values in any of the predictors
    formula.onlypreds <- formula[-2]
    
    X <- model.matrix(formula.onlypreds,
      model.frame(formula.onlypreds, datanomissboth,
        na.action = stats::na.omit))
    
    ## divide data set into sampled sites and unsampled sites based
    ## on whether the response variable has a number (for sampled sites)
    ## or NA (for unsampled sites)
    
    ind.sa <- !is.na(yvar)
    ind.un <- is.na(yvar)
    data.sa <- datanomissboth[ind.sa, ]
    data.un <- datanomissboth[ind.un, ]
    
    m.un <- stats::model.frame(formula, data.un, na.action =
        stats::na.pass)
    
    
    Xu <- model.matrix(formula.onlypreds,
      model.frame(formula.onlypreds, data.un,
        na.action = stats::na.omit))
    
    
    ## sampled response values and design matrix
    m.sa <- stats::model.frame(formula, data.sa, na.action =
        stats::na.omit)
    w.sa <- stats::model.response(m.sa)
    Xs <- stats::model.matrix(formula, m.sa)
    w.density <- w.sa
    n <- nrow(Xs)
    
    prednames <- colnames(Xs)
    
    ## x and y coordinates for sampled and unsampled sites
    x.sa <- xcoordsUTM[ind.sa]
    y.sa <- ycoordsUTM[ind.sa]
    x.un <- xcoordsUTM[ind.un]
    y.un <- ycoordsUTM[ind.un]
    
    ## number of sites that were sampled
    n.sa <- nrow(Xs)
    ## number of sites that were not sampled
    n.un <- nrow(Xu)
    
    
    ## estimate the spatial parameters, the covariance matrix, and
    ## the inverse of the covariance matrix
    
    spat.est <- estcovparm(response = density,
      designmatrix = as.matrix(X),
      xcoordsvec = xcoordsUTM,
      ycoordsvec = ycoordsUTM, CorModel = CorModel,
      estmethod = "ML",
      covestimates = c(NA, NA, NA),
      pivec = piest,
      Vnn = Vnn)
    
    parms.est <- spat.est$parms.est
    Sigma <- spat.est$Sigma
    min2loglik <- spat.est$min2loglik
    Sigma.ss <- Sigma[ind.sa, ind.sa]
    
    nugget.effect <- parms.est[1]; parsil.effect <- parms.est[2]
    range.effect <- parms.est[3]
    
    
    Sigma.ssi <- NA
    
    ## the generalized least squares regression coefficient estimates
    
    betahat <- spat.est$b.hat
    
    ## estimator for the mean vector
    muhats <- as.matrix(Xs) %*% betahat
    muhatu <- as.matrix(Xu) %*% betahat
    
    resids <- w.sa - muhats * piest[ind.sa]
    
    muhat <- rep(NA, nrow(datanomissboth))
    muhat[ind.sa == TRUE] <- muhats
    muhat[ind.sa == FALSE] <- muhatu
    
    ## construct matrices necessary for FPBK estimator using
    ## detection.
    
    Sigma.su <- Sigma[ind.sa, ]
    D <- Sigma
    R <- piest[ind.sa] * Sigma.su
    
    C <- matrixcalc::hadamard.prod(as.matrix(piest[ind.sa]) %*%
        t(as.matrix(piest[ind.sa])),
      Sigma.ss) +
      matrixcalc::hadamard.prod(as.matrix(muhats) %*%
          t(as.matrix(muhats)), Vnn[ind.sa, ind.sa]) +
      matrixcalc::hadamard.prod(Sigma.ss, Vnn[ind.sa, ind.sa]) +
      diag(as.vector(muhats) *
          (as.vector(piest[ind.sa])) *
          (1 - as.vector(piest[ind.sa])),
        nrow = nrow(Sigma.ss))

    Cinv <- solve(C)

    Xnstar <- piest[ind.sa] * Xs
    
    ## returns a list with the following components:
    ## 1.) A vector of the estimated regression coefficients
    ## 2.) A vector of the estimated covariance parameters
    ## 3.) A list with the information needed by FPBKpred
    
    covparms <- as.vector(c(nugget.effect, parsil.effect,
      range.effect))
    betahatest <- as.vector(betahat)
    covest <- spat.est$covb
    
    
    names(covparms) <- c("Nugget", "Partial Sill", "Range")
    
    FPBKpredobj <- list(formula, datanomissboth, xcoordsUTM,
      ycoordsUTM,
      estmethod, CorModel, Sigma, Sigma.ssi, C,
      Cinv, Xnstar, R, piest)
    names(FPBKpredobj) <- c("formula", "data", "xcoordsUTM",
      "ycoordsUTM",
      "estmethod","correlationmod", "covmat", "covmatsampi", "C",
      "Cssi", "Xnstar", "R", "piest")
    obj <- list(covparms, betahatest, covest, min2loglik, prednames,
      n, CorModel, resids, Xs, w.sa, detind, FPBKpredobj)
    
    names(obj) <- c("SpatialParmEsts", "CoefficientEsts",
      "BetaCov", "minus2loglike", "PredictorNames", "SampSize",
      "CovarianceMod",
      "resids", "DesignMat", "Density", "DetectionInd",
      "FPBKpredobj")
    
  }
  class(obj) <- "slmfit"
  return(obj)
  
  }

##counts <- 1:10
##pred1 <- runif(10, 0, 1)
##pred2 <- runif(10, 0, 2)
##data <- data.frame(cbind(counts, pred1, pred2))
##formula <- counts ~ pred1 + pred2
#

#  summary(object = slm_info)
# print(x = summary(object = slm_info))
# print(slm_info)
# summary(slm_info)


# missingind <- sample(1:length(simdata$X5), size = 40)
#
# X5mis <- simdata$X5
# X5mis[missingind] <- NA
# simdata$X5mis <- X5mis
# zmis <- simdata$Z
# zmis[sample(1:400, size = 71)] <- NA
# simdata$zmis <- zmis
# slm_info <- slmfit(zmis ~ X1 + X2 + X5mis, data = simdata,
#  xcoordcol = "x", ycoordcol = "y",  coordtype = "UTM",
#    estmethod = "ML")


#designmatrixsa <- with(exampledataset, model.matrix(counts ~ pred1 + pred2))
##ind.sa <- is.na(exampledataset$counts) == FALSE
##response <- exampledataset$counts
# solve(t(Xtest) %*% Xtest) %*% t(Xtest) %*% matrix(exampledataset$counts[is.na(exampledataset$counts) == FALSE])

##print.summary.slmfit(x = summary.slmfit(object = slm_info))
##summary(slm_info)
##print(slm_info)

##pred_info <- predict(object = slm_info, FPBKcol = NULL)

##pred_info

##FPBKoutput(pred_info = pred_info, get_variogram = TRUE,
##  get_sampdetails = TRUE,
##  get_krigmap = FALSE, get_report = TRUE, conf_level = c(0.80, 0.90,
##    0.95))