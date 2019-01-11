#' Estimate Covariance Parameters
#'
#' Used to estimate spatial covariance parameters for a few different spatial models.
#' Estimated parameters can then be used in `FPBKpred` to predict unobserved values.
#'
#' The parameters used in the function can be inherited from the FPBKpred function
#'
#' @param formula is  an R linear model formula specifying the counts as the
#' response variable as well as covariates for predicting counts of animals or plants
#' on the unsampled sites.
#' @param data is the data set with the response column of counts, the covariates to
#' be used for the block kriging, and the spatial coordinates for all of the sites.
#' @param xcoordcol is the name of the column in the data frame with x coordinates
#' @param ycoordcol is the name of the column in the data frame with y coordinates
#' @param CorModel is the covariance structure. By default, `covstruct` is
#' Exponential but other options include the Matern, Spherical, and Gaussian.
#' @return a list with a vector of estimated covariance parameters, the fitted
#' covariance matrix, and the inverse of the fitted covariance matrix.
#' @export estcovparm

estcovparm <- function(formula, data, xcoordcol, ycoordcol,
  CorModel = "Exponential") {

  ## only estimate parameters using sampled sites only
  fullmf <- stats::model.frame(formula, na.action = 
      stats::na.pass)
  yvar <- stats::model.response(fullmf, "numeric")
  
  ind.sa <- !is.na(yvar)
  data.sa <- data[ind.sa, ]


  names.theta <- c("nugget", "parsil", "range")

  ## eventually will expand this section to include other covariance types
  if(CorModel == "Cauchy" || CorModel == "BesselK"){
    names.theta <- c(names.theta, "extrap")
  }

  nparm <- length(names.theta)

  n <- length(yvar[ind.sa])
  p <- length(attr(stats::terms(formula), "variables")) - 2


  ## distance matrix for all of the sites
  distmatall <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  distmatall[lower.tri(distmatall)] <- stats::dist(as.matrix(cbind(data$xcoordsUTM, data$ycoordsUTM)))
  distmatall <- distmatall + t(distmatall)

  ## constructing the distance matrix between sampled sites only
  sampdistmat <- matrix(0, n, n)
  sampdistmat[lower.tri(sampdistmat)] <- stats::dist(as.matrix(cbind(data$xcoordsUTM[ind.sa], data$ycoordsUTM[ind.sa])))
  distmat <- sampdistmat + t(sampdistmat)

  ## perform a grid search on log scale to find an appropriate
  ## starting value for optim. The grid search covers a variety
  ## of nugget to partial sill ratios as well as a few different
  ## range parameters spanning the maximum distance of the distance matrix

  possible.nugget <- c(stats::var(yvar[ind.sa]),
    stats::var(yvar[ind.sa]) / 2,
    stats::var(yvar[ind.sa]) / 4)
  possible.theta1 <- log(possible.nugget)
  possible.parsil <- c(stats::var(yvar[ind.sa]),
    stats::var(yvar[ind.sa]) / 2,
    stats::var(yvar[ind.sa]) / 4)
  possible.theta2 <- log(possible.parsil)
  possible.range <- c(max(distmat), max(distmat) / 2, max(distmat) / 4)
  possible.theta3 <- log(possible.range)

  theta <- expand.grid(possible.theta1, possible.theta2, possible.theta3)


  ## construct the design matrix based on the formula input
  XDesign <- stats::model.matrix(formula)

  m2loglik <- rep(NA, nrow(theta))

if (CorModel == "Exponential") {
    for (i in 1:nrow(theta)) {
      m2loglik[i] <- m2LL.FPBK.nodet.exp(theta = theta[i, ], zcol = yvar[ind.sa],
        XDesign = XDesign,
        xcoord = data.sa$xcoordsUTM, ycoord = data.sa$ycoordsUTM)
    }

    max.lik.obs <- which(m2loglik == min(m2loglik))

    ## optimize using Nelder-Mead
    parmest <- optim(theta[max.lik.obs, ], m2LL.FPBK.nodet.exp,
      zcol = yvar[ind.sa],
      XDesign = XDesign,
      xcoord = data.sa$xcoordsUTM, ycoord = data.sa$ycoordsUTM,
      method = "Nelder-Mead")

    ## extract the covariance parameter estimates. When we deal with covariance
    ## functions with more than 3 parameters, this section will need to be modified
    min2loglik <- parmest$value

    nugget.effect <- exp(parmest$par[1])
    parsil.effect <- exp(parmest$par[2])
    range.effect <- exp(parmest$par[3])
    parms.est <- c(nugget.effect, parsil.effect, range.effect)

    Sigma <- diag(nugget.effect,
      nrow = nrow(distmatall)) +
      parsil.effect * exp(-distmatall / range.effect)

    Sigmai <- solve(Sigma)

} else if (CorModel == "Spherical") {

  for (i in 1:nrow(theta)) {
    m2loglik[i] <- m2LL.FPBK.nodet.sph(theta = theta[i, ], zcol = yvar[ind.sa],
      XDesign = XDesign,
      xcoord = data.sa$xcoordsUTM, ycoord = data.sa$ycoordsUTM)
  }

  max.lik.obs <- which(m2loglik == min(m2loglik))

  ## optimize using Nelder-Mead
  parmest <- optim(theta[max.lik.obs, ], m2LL.FPBK.nodet.sph,
    zcol = data.sa$counts,
    XDesign = XDesign,
    xcoord = data.sa$xcoordsUTM, ycoord = data.sa$ycoordsUTM,
    method = "Nelder-Mead")

  min2loglik <- parmest$value

  nugget.effect <- exp(parmest$par[1])
  parsil.effect <- exp(parmest$par[2])
  range.effect <- exp(parmest$par[3])
  parms.est <- c(nugget.effect, parsil.effect, range.effect)

  cormatSpher <- 1 - (3 / 2) * (distmat / range.effect) +
    (1 / 2) * (distmat / range.effect) ^ 3
  cormatSpher[distmat > range.effect] <- 0
  Sigma <- diag(nugget.effect, nrow = nrow(distmat)) +
    parsil.effect * cormatSpher

  Sigmai <- solve(Sigma)
}

  return(list(parms.est, Sigma, Sigmai))
}

counts <- c(1, NA, NA, NA, 3, 1:13, 21, 30)
pred1 <- runif(20, 0, 1); pred2 <- rnorm(20, 0, 1)
xcoords <- runif(20, 0, 1); ycoords <- runif(20, 0, 1)
dummyvar <- runif(20, 0, 1)
xcoordcol <- "xcoords"; ycoordcol <- "ycoords"
CorModel = "Exponential"

data <- as.data.frame(cbind(counts, pred1, pred2, xcoords, ycoords, dummyvar))

formula <- counts ~ pred1 + pred2

##estcovparm(formula = formula, data = data, xcoordcol = xcoordcol,
##  ycoordcol = ycoordcol, CorModel = "Exponential")[[1]]
