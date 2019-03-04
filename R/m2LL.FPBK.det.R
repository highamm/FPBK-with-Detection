#' Covariance Parameter Estimation Function.
#'
#' The primary purpose of \code{m2LL.FPBK.det} is to estimate the spatial
#' covariance parameters using Maximum Liklihood when detection
#' is included.
#'
#' @param theta is the parameter vector of (nugget, partialsill, range)
#' @param zcol is the response vector of counts
#' @param XDesign is the design matrix containing the covariates used to predict animal or plant abundance (including a column of 1's for the intercept).
#' @param xcoord is a vector of the x spatial coordinates (in UTM)
#' @param ycoord is a vector of the y spatial coordinates (in UTM)
#' @param CorModel is the geostatistical spatial correlation model to be used. See the \code{corModels} documentation for possible models to use.

#' @return A numeric output of minus 2 times the restricted log likelihood to be minimized by `optim` to obtain spatial parameter estimates.
#' @importFrom stats optim
#' @importFrom stats glm
#' @importFrom stats rbinom
#' @export m2LL.FPBK.det

## split into different functions for different covariance matrix structures

m2LL.FPBK.det <- function(theta, zcol, XDesign, xcoord, ycoord,
  CorModel, pivec, Vnn) {
  
  n <- length(zcol)
  p <- length(XDesign[1, ])
  nugget <- as.numeric(exp(theta[1]))
  parsil <- as.numeric(exp(theta[2]))
  range <- as.numeric(exp(theta[3]))
  mu <- as.numeric(exp(theta[4]))
  ##nugget.z <- as.numeric(exp(theta[5]))
  DM <- matrix(0, n, n)
  DM[lower.tri(DM)] <- stats::dist(as.matrix(cbind(xcoord, ycoord)))
  Dismat <- DM + t(DM)
  Sigmat <- parsil * exp(-Dismat / range)
  Cmat <- diag(as.vector(mu) *
      ##(as.vector(pivec) - diag(Vnn) - as.vector(pivec^2)),
      (as.vector(pivec) * (1 - as.vector(pivec))),
    nrow = nrow(Sigmat)) + 
    (pivec %*% t(pivec)) * (diag(nugget, nrow = nrow(Sigmat)) + Sigmat) +
    (as.matrix(rep(mu, n)) %*% t(as.matrix(rep(mu, n)))) * Vnn  +
    (diag(nugget, nrow = nrow(Sigmat)) + Sigmat) * Vnn
  
  Ci <- solve(Cmat, tol = 1e-27)
  minusloglik <- (1 / 2) * determinant(Cmat, logarithm = TRUE)$modulus +
    (1 / 2) * (t(as.matrix(zcol) -
        as.matrix(rep(mu, n) * as.vector(pivec)))) %*% Ci %*%
    (as.matrix(zcol) - as.matrix(rep(mu, n) * as.vector(pivec)))
  ##covbi <- t(XDesign) %*% Ci %*% XDesign
  ##covb <- solve(covbi, tol = 1e-21)
  ##b.hat <- covb %*% t(XDesign) %*% Ci %*% zcol
  ##r <- zcol - XDesign %*% b.hat
  ##t(r) %*% Ci %*% r + sum(log(svd(Cmat)$d)) + sum(log(svd(covbi)$d))
  return(as.numeric(minusloglik))
}
