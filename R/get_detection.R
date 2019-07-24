#' Obtains estimates for detection parameters through logistic regression.
#'
#' Estimates logistic regression coefficients from detection data.
#'
#' @param formula is an R linear model formula specifying binary detection as the
#' response variable as well as covariates associated with detection.
#' @param data is the data set with the response column of detection and the covariates to
#' be used in the logistic regression, and the spatial coordinates for all of the sites.
#' @param varmethod is the method to obtain the variance matrix for the detection probabilities, either "Delta" or "Bootstrap".
#' @return a list with \itemize{
#'   \item Fitted detection coefficients
#'   \item The asymptotic variance-covariance matrix of the 
#'   fitted regression coefficients
#'   \item the formula used to fit the logistic regression model
#'   \item A list of bootstrapped regression coefficients used to 
#'   obtain a covariance matrix for detection estimates
#'   \item The method used to get the variance covariance matrix
#'   for detection probabilities.
#' }
#' @examples
#' df <- data.frame(y = rbinom(30, 1, 0.7), x = rnorm(30, 4, 0.4))
#' get_detection(formula = y ~ x, data = df)
#' @import stats
#' @export get_detection

get_detection <- function(formula, data, 
  varmethod = "Bootstrap") {
  
  Xall <- model.matrix(formula, model.frame(formula, data,
    na.action = stats::na.pass))
  
  missingind <- base::apply(is.na(Xall), MARGIN = 1, FUN = sum)
  nmissing <- sum(missingind >= 1)
  
  datanomiss <- as.data.frame(data[missingind == 0, ])
  colnames(datanomiss) <- colnames(data)
  
  if (nmissing >= 1) {
    warning(paste("There were", nmissing, "sites with predictors with missing values. These will be removed from the data set and further analysis will be completed without these observations."))
  }
  
  fullmf <- stats::model.frame(formula, na.action =
      stats::na.pass, data = datanomiss)
  detection_resp <- stats::model.response(fullmf, "numeric")
  
  if (length(detection_resp) < 30) {
    warnings("It is recommended that the size of the detection data set
      is larger than 30 for the logistic regression, particularly 
      if detection is small")
  } 
  
  X <- model.matrix(formula,
    model.frame(formula, datanomiss,
      na.action = stats::na.omit))
  
  mod <- stats::glm(detection_resp ~ X - 1, family = "binomial")
  coefs <- base::summary(mod)$coef[, 1]
  
  F <- base::summary(mod)$cov.unscaled
  fitted.values <- mod$fitted.values
  
  if (varmethod == "Bootstrap") {
  niter <- 1400
  boot.ind.gen <- function() {
    boot.ind <- sample(1:nrow(X), size = nrow(X), replace = TRUE)
    return(boot.ind)
  }
  
  ind.res <- replicate(niter, boot.ind.gen())
  boot.coefs <- matrix(NA, nrow = niter, 
    ncol = length(coefs))
  
  for (kk in 1:ncol(ind.res)) {
    
    boot.mat <- cbind(detection_resp, X)[ind.res[ ,kk], ]
    
    if (mean(boot.mat[ ,1]) == 1) {
      boot.coefs[kk, ] <- c(10, rep(0, ncol(X) - 1))
      warning("In bootstrapping, one iteration had all responses
        equal to 1 so detection coefficients were set to an 
        intercept of 10 and predictor coefficients of 0 for
        this iteration.")
    } else if (mean(boot.mat[ ,1]) == 0) {
      boot.coefs[kk, ] <- c(-10, rep(0, ncol(X) - 1))
      warning("In bootstrapping, one iteration had all responses
        equal to 0 so detection coefficients were set to an 
        intercept of -10 and predictor coefficients of 0 for
        this iteration.")
    } else{
    
    ## fit the logistic regression with the resampled rows
      if (ncol(boot.mat) > 2){
    boot.mod <- base::suppressWarnings(stats::glm(boot.mat[ ,1] ~
        boot.mat[ ,-c(1, 2)],
      family = "binomial"))
      } else {
        boot.mod <- base::suppressWarnings(stats::glm(boot.mat[ ,1] ~
            1,
          family = "binomial"))
      }
    boot.coefs[kk, ] <- base::suppressWarnings(base::summary(boot.mod)$coef[ ,1])
    
    }
    
  }
  
  } else if (varmethod == "Delta") {
    boot.coefs <- NULL
  } else {
    stop("varmethod must be either 'Delta' or 'Bootstrap' to obtain
      covariance matrix for detection probabilities")
  }
  
  
  
  
  obj <- list(coefs = coefs, covmat = F, formula = formula,
    boot.coefs = boot.coefs, varmethod = varmethod)
  
  class(obj) <- "get_detection"
  return(obj)
  
}




