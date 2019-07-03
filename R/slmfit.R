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
#' @param estmethod is either the default \code{"REML"} for restricted
#' maximum likelihood to estimate the covariance parameters and
#' regression coefficients or \code{"ML"} to estimate the covariance
#' parameters and regression coefficients.
#' @param covestimates is an optional vector of covariance parameter estimates (nugget, partial sill, range). If these are given and \code{estmethod = "None"}, the the provided vector are treated as the estimators to create the covariance structure.
#' @param detectionobj is a fitted model obj from \code{get_detection}. The default is for this object to be \code{NULL}, resulting in
#' spatial prediction that assumes perfect detection.
#' @param areacol is the name of the column with the areas of the sites. By default, we assume that all sites have equal area, in which
#' case a vector of 1's is used as the areas.
#' @return a list with \itemize{
#'   \item the spatial covariance estimates
#'   \item the regression coefficient estimates
#'   \item the estimated covariance matrix
#'   \item minus two times the log-likelihood
#'   \item the names of the predictors used as covariates
#'   \item the sample size
#'   \item the name of the correlation model used
#'   \item a vector of residuals
#'   \item the design matrix for the sampled sites
#'   \item a vector of the density of the observed counts
#'   \item an indicator for whether the user input detection data
#'   \item a list with information for \code{predict} containing \enumerate{
#'        \item formula
#'        \item data, the data set without any missing values
#'        \item a vector of x-coordinates (in UTM)
#'        \item a vector of y-coordinates (in UTM)
#'        \item whether REML or ML was used
#'        \item the correlation model used
#'        \item the covariance matrix for all of the sites
#'        \item Inverted covariance matrix on the sampled sites
#'        \item Other items used in \code{predict}
#'        }
#' }
#' @import stats
#' @export slmfit

slmfit <- function(formula, data, xcoordcol, ycoordcol,
  CorModel = "Exponential",
  estmethod = "REML",
  covestimates = c(NA, NA, NA),
  detectionobj = NULL,
  areacol = NULL) {
  
  ## make sure that detectionobj comes from get_detection if included
  
  if (is.null(detectionobj) == FALSE) {
    if (class(detectionobj) != "get_detection") {
      stop("detectionobj must be of class 'get_detection', generated
        from the get_detection function.")
    }
  }
  
  ## make sure estmethod is either REML, ML, or None
  
  if (estmethod != "REML" & estmethod != "ML" &
      estmethod != "None") {
    stop("estmethod must be either 'REML' for restricted maximum
      likelihood, 'ML' for maximum likelihood, or 'None' with
      covariance parameters specified in the covestimates
      argument.")
  }
  ## make sure CorModel is one of the options we have set-up
  if (CorModel != "Exponential" & CorModel != "Spherical" &
      CorModel != "Gaussian") {
    stop("'CorModel' must be either 'Exponential', 'Spherical', or
      'Gaussian'")
  }
  
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
  
  if (estmethod == "None" & sum(covestimates < 0) > 0) {
    stop("'covestimates' must be a vector of positive values with
      the (nugget, partial sill, range) specified.")
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
  
  if (sum(names(data) == xcoordcol) == 0 |
      sum(names(data) == ycoordcol) == 0) {
    stop("xcoordcol and ycoordcol must be the names of the columns
      in the data set (in quotes) that specify the x and y coordinates.")
  }
  
  ## convert all character predictor variables into factors,
  ## with a warning message.
  # datapredsonly <- data.frame(data[ ,attr(terms(formula), "term.labels")])
  # colnames(datapredsonly) <- attr(terms(formula), "term.labels")
  # predictormatch <- match(names(data), names(datapredsonly))
  # 
  # if (ncol(datapredsonly) >= 1) {
  # 
  # if (sum(sapply(datapredsonly, is.character)) > 0) {
  #   warning("At least one predictor variable is a character, which has been converted into a factor.")
  # }
  # 
  # data[ ,sapply(data, is.character) & is.na(predictormatch) == FALSE] <- factor(data[ ,which(sapply(data, is.character))])
  # 
  # 
  # ## check to make sure number of factor levels is somewhat small.
  # ## If not, return a warning.
  # if (max(sapply(datapredsonly[ ,sapply(datapredsonly, is.factor)], nlevels)) > 20) {
  #   warning("At least one predictor variable has more than 20 factor levels.")
  # }
  # 
  # }
  Xall <- model.matrix(formula, model.frame(formula, data,
    na.action = stats::na.pass))
  
  missingind <- base::apply(is.na(Xall), MARGIN = 1, FUN = sum)
  nmissing <- sum(missingind >= 1)
  
  datanomiss <- data[missingind == 0, ]
  
  ## give a warning if some of the predictors have missing values.
  if (nmissing >= 1) {
    warning(paste("There were", nmissing, "sites with predictors with missing values. These will be removed from the data set and further analysis will be completed without these observations."))
  }
  
  ## ASSUME that coordinates are in TM

    xcoordsUTM <- datanomiss[ ,xcoordcol]
    ycoordsUTM <- datanomiss[ ,ycoordcol]
  
  
  detind <- is.null(detectionobj)
  
  if (is.null(areacol) == FALSE) {
    if (is.numeric(datanomiss[ ,areacol]) == FALSE |
        sum(is.na(datanomiss[ ,areacol])) > 0) {
      stop("'areacol' must specify the name of the column in the data set with the areas for each site. This column must be numeric 
        without any missing values.")
    }
  }
  if (is.null(detectionobj) == TRUE) {
  ## create the design matrix for unsampled sites, for all of the sites, and for the sampled sites, respectively.
  
    
    
    if (is.null(areacol) == TRUE) {
      areavar <- rep(1, nrow(datanomiss))
    } else {
      areavar <- datanomiss[ ,areacol]
    }
    
  fullmf <- stats::model.frame(formula, na.action =
      stats::na.pass, data = datanomiss)
  yvar <- stats::model.response(fullmf, "numeric")
  density <- yvar / areavar
  
  if (is.numeric(yvar) == FALSE) {
    stop("Check to make sure response variable is numeric, not a factor or character.")
  }
  
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
  z.density <- z.sa / areavar[ind.sa]
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
  
  resids <- z.density - muhats
  
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
    estmethod, CorModel, Sigma, Sigma.ssi, areavar)
  names(FPBKpredobj) <- c("formula", "data", "xcoordsUTM",
    "ycoordsUTM",
    "estmethod","correlationmod", "covmat", "covmatsampi",
    "areavar")
  obj <- list(covparms, betahatest, covest, min2loglik, prednames,
    n, CorModel, resids, Xs, z.density, detind, FPBKpredobj)
  
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
    piboot <- NULL
    
    } else if (detectionobj$varmethod == "Bootstrap") {
      
    boot.coefs <- detectionobj$boot.coefs
    logitboot <- Xndet %*% t(boot.coefs)
    piboot <- exp(logitboot) / (1 + exp(logitboot))
    Vnn <- cov(t(piboot), use = "complete.obs")
    
    }
    
    if (is.null(areacol) == TRUE) {
      areavar <- rep(1, nrow(datanomissboth))
    } else {
      areavar <- datanomissboth[ ,areacol]
    }
    
    fullmf <- stats::model.frame(formula, na.action =
        stats::na.pass, data = datanomissboth)
    yvar <- stats::model.response(fullmf, "numeric")
    density <- yvar / areavar
    
    if (is.numeric(yvar) == FALSE) {
      stop("Check to make sure response variable is numeric, not a factor or character.")
    }
    
    ## remove any rows with missing values in any of the predictors
    formula.onlypreds <- formula[-2]
    
    X <- model.matrix(formula.onlypreds,
      model.frame(formula.onlypreds, datanomissboth,
        na.action = stats::na.omit))
    
    ## divide data set into sampled sites and unsampled sites based
    ## on whether the response variable has a number (for sampled sites)
    ## or NA (for unsampled sites)
    
    ind.sa <- !is.na(density)
    ind.un <- is.na(density)
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
    w.density <- w.sa / areavar[ind.sa]
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
    
    resids <- w.density - muhats * piest[ind.sa]
    
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

    Cinv <- solve(C, tol = 1e-23)

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
      Cinv, Xnstar, R, piest, areavar)
    names(FPBKpredobj) <- c("formula", "data", "xcoordsUTM",
      "ycoordsUTM",
      "estmethod","correlationmod", "covmat", "covmatsampi", "C",
      "Cssi", "Xnstar", "R", "piest", "areavar")
    obj <- list(covparms, betahatest, covest, min2loglik, prednames,
      n, CorModel, resids, Xs, w.density, detind, FPBKpredobj,
      piboot)
    
    names(obj) <- c("SpatialParmEsts", "CoefficientEsts",
      "BetaCov", "minus2loglike", "PredictorNames", "SampSize",
      "CovarianceMod",
      "resids", "DesignMat", "Density", "DetectionInd",
      "FPBKpredobj", "piboot")
    
  }
  class(obj) <- "slmfit"
  return(obj)
  
  }

