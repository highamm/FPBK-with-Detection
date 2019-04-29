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
#' @param pivec is a vector of estimated detection probabilities
#' on each of the sites
#' @param Vnn is the covariance matrix for the estimated detection
#' probabilities
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
  mu <- XDesign %*% as.matrix(theta[4:length(theta)])
  
  ##nugget.z <- as.numeric(exp(theta[5]))
  DM <- matrix(0, n, n)
  DM[lower.tri(DM)] <- stats::dist(as.matrix(cbind(xcoord, ycoord)))
  Dismat <- DM + t(DM)
  
  if (CorModel == "Exponential") {
    Sigmat <- parsil * exp(-Dismat / range)
  } else if (CorModel == "Gaussian") {
    Sigmat <- parsil * exp(-Dismat ^ 2 / range)
    } else if (CorModel == "Spherical") {
      ## CHANGE THIS
      Sigmat <- (1 - 1.5 * (Dismat / range) +
          0.5 * (Dismat / range) ^ 3)
      Sigmat[Dismat / range > 1] <- 0  
      Sigmat <- parsil * Sigmat
      }
  
  Cmat <- diag(as.vector(mu) *
      ##(as.vector(pivec) - diag(Vnn) - as.vector(pivec^2)),
      (as.vector(pivec) * (1 - as.vector(pivec))),
    nrow = nrow(Sigmat)) + 
    (pivec %*% t(pivec)) * (diag(nugget, nrow = nrow(Sigmat)) + Sigmat) +
    (as.matrix(mu) %*% t(as.matrix(mu))) * Vnn  +
    (diag(nugget, nrow = nrow(Sigmat)) + Sigmat) * Vnn
  
  nug_prop <- nugget / (nugget + parsil)
  if (nug_prop < 0.0001) {
    Cmat <- Cmat + diag(1e-6, nrow = nrow(Cmat))
  }
  
  Ci <- solve(Cmat, tol = 1e-27)
  minusloglik <- (1 / 2) * determinant(Cmat, logarithm = TRUE)$modulus +
    (1 / 2) * (t(as.matrix(zcol) -
        as.matrix(mu * as.vector(pivec)))) %*% Ci %*%
    (as.matrix(zcol) - as.matrix(mu * as.vector(pivec)))
  ##covbi <- t(XDesign) %*% Ci %*% XDesign
  ##covb <- solve(covbi, tol = 1e-21)
  ##b.hat <- covb %*% t(XDesign) %*% Ci %*% zcol
  ##r <- zcol - XDesign %*% b.hat
  ##t(r) %*% Ci %*% r + sum(log(svd(Cmat)$d)) + sum(log(svd(covbi)$d))
  return(as.numeric(minusloglik))
}
