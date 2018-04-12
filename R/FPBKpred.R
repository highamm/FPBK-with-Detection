#' Perform Finite Population Block Kriging
#'
#' Takes sample data and uses FPBK to predict the counts on the unsampled sites.
#' The column with the counts should have numeric values for the observed counts
#' on the sampled sites and `NA` for any site that was not sampled.
#'
#' @param formula is an R linear model formula specifying the counts as the
#' response variable as well as covariates for predicting counts of animals or plants
#' on the unsampled sites.
#' @param data is the data set with the response column of counts, the covariates to
#' be used for the block kriging, and the spatial coordinates for all of the sites.
#' @param xcoordcol is the name of the column in the data frame with x coordinates
#' @param ycoordcol is the name of the column in the data frame with y coordinates
#' @param covstruct is the covariance structure. By default, `covstruct` is
#' Exponential but other options include the Matern, Spherical, and Gaussian.
#' @param FPBK.col is a vector in the data set that contains the weights for
#' prediction. The default setting predicts the population total
#' @return a list with the estimated population total, the estimated prediction
#' variance, and the vector of predicted counts for all of the sites.
#' @export FPBKpred



FPBKpred <- function(formula, data, xcoordcol, ycoordcol,
  covstruct = "Exponential", FPBK.col = NULL) {

  ## if FPBK.col is left out, we are predicting the population total.
  ## Otherwise, FPBK.col is the name of the column in the data set
  ## with the weights for the sites that we are predicting (eg. a vector
  ## of 1's and 0's for predicting the total of the sites marked with 1's)

  if (is.null(FPBK.col) == TRUE) {
    predwts <- rep(1, nrow(data))
  } else if (is.character(FPBK.col) == TRUE) {
    predwts <- data[ ,FPBK.col]
  } else{
    stop("FPBK.col must be a character specifying the name of the column of
      prediction weights in the data set")
  }

  ## divide data set into sampled sites and unsampled sites based
  ## on whether the response variable has a number (for sampled sites)
  ## or NA (for unsampled sites)
  response.col <- as.character(attr(stats::terms(formula, data = data), "variables"))[2]
  ind.sa <- !is.na(data[, response.col])
  ind.un <- is.na(data[, response.col])
  data.sa <- data[ind.sa, ]
  data.un <- data[ind.un, ]

  ## change factor count responses or character count responses to numeric
  if(is.factor(data.sa[ ,response.col])) {
    data.sa[ ,response.col] <- as.numeric(as.character(data.sa[ ,response.col]))
  }
  if(is.character(data.sa[ ,response.col])) {
    data.sa[ ,response.col] <- as.numeric(data.sa[ ,response.col])
  }

  ## display some warnings if the user, for example tries to input the
  ## vector of xcoordinates as the input instead of the name of the column
  if(is.character(xcoordcol) == FALSE) {
    stop("xcoords must be a character giving the name of the column in the data set
      with the x coordinates")
  }
  if(is.character(ycoordcol) == FALSE) {
    stop("ycoords must be a character giving the name of the column in the data set
      with the y coordinates")
  }

  ## put a 1 here for each value of the response in order to create the model matrix
  data.un[ ,response.col] <- 1

  ## create the design matrix for unsampled sites, for all of the sites,
  ## and for the sampled sites, respectively.
  m.un <- stats::model.frame(formula, data.un)
  Xu <- stats::model.matrix(formula, m.un)

  X <- stats::model.matrix(formula, stats::model.frame(formula,
    data, na.action = "na.pass"))

  ## sampled response values
  z.sa <- stats::model.frame(formula, data.sa)
  z.sa <- matrix(data.sa[ ,names(z.sa)[1]], ncol = 1)
  m.sa <- stats::model.frame(formula, data.sa)
  Xs <- stats::model.matrix(formula, m.sa)

  ## x and y coordinates for sampled and unsampled sites
  x.sa <- data.sa[ ,xcoordcol]
  y.sa <- data.sa[ ,ycoordcol]
  x.un <- data.un[ ,xcoordcol]
  y.un <- data.un[ ,ycoordcol]

  ## number of sites that were sampled
  n.sa <- sum(ind.sa)
  ## number of sites that were not sampled
  n.un <- sum(ind.un)

  B <- predwts
  Bs <- B[ind.sa]
  Bu <- B[ind.un]

  ## estimate the spatial parameters, the covariance matrix, and
  ## the inverse of the covariance matrix
  spat.est <- estcovparm(formula = formula, data = data, xcoordcol = xcoordcol,
    ycoordcol, CorModel = covstruct)
  parms.est <- spat.est[[1]]
  Sigma <- spat.est[[2]]
  Sigmai <- spat.est[[3]]
  nugget.effect <- parms.est[1]; parsil.effect <- parms.est[2]
  range.effect <- parms.est[3]

  ## used in the Kriging formulas
  Sigma.us <- Sigma[ind.un, ind.sa]
  Sigma.su <- t(Sigma.us)
  Sigma.ss <- Sigma[ind.sa, ind.sa]
  Sigma.ssi <- solve(Sigma.ss)

  ## the generalized least squares regression coefficient estimates
  betahat <- solve((t(Xs) %*% Sigma.ssi %*% Xs)) %*%
     (t(Xs) %*% Sigma.ssi %*% as.matrix(data.sa[ ,response.col]))

  ## estimator for the mean vector
  muhats <- Xs %*% betahat
  muhatu <- Xu %*% betahat

  ## matrices used in the kriging equations
  Cmat <- Sigma.ss %*% as.matrix(Bs) + Sigma.su %*% as.matrix(Bu)
  Dmat <- t(X) %*% matrix(B) - t(Xs) %*% Sigma.ssi %*% Cmat
  Vmat <- solve(t(Xs) %*% Sigma.ssi %*% Xs)

  ## the predicted values for the sites that were not sampled
  zhatu <- Sigma.us %*% Sigma.ssi %*% (z.sa -
      muhats) + muhatu

  ## creating a column in the outgoing data set for predicted counts as
  ## well as a column indicating whether or not the observation was sampled
  ## or predicted
  data$preds <- data[ ,response.col]
  data$preds[is.na(data$preds) == TRUE] <- zhatu
  data$sampind <- 1
  data$sampind[is.na(data[ ,response.col]) == TRUE] <- 0

  ## the FPBK predictor
  FPBKpredictor <- t(Bs) %*% z.sa + t(Bu) %*% zhatu

  ## the prediction variance for the FPBK predictor
  pred.var.obs <- t(as.matrix(B)) %*% Sigma %*%
    as.matrix(B) -
    t(Cmat) %*% Sigma.ssi %*% Cmat +
    t(Dmat) %*% Vmat %*% Dmat

  return(list(FPBKpredictor, pred.var.obs,
    as.matrix(cbind(data[ ,xcoordcol], data[, ycoordcol], data$preds, data$sampind)),
    as.vector(c(nugget.effect, parsil.effect, range.effect))))
}


counts <- c(1, NA, NA, NA, 3, 1:35)
pred1 <- runif(40, 0, 1); pred2 <- rnorm(40, 0, 1)
xcoords <- runif(40, 0, 1); ycoords <- runif(40, 0, 1)
dummyvar <- runif(40, 0, 1)
df <- as.data.frame(cbind(counts, pred1, pred2, xcoords, ycoords, dummyvar))
data <- df
FPBK.col <- NULL
xcoordcol <- "xcoords"; ycoordcol <- "ycoords"
formula <- counts ~ pred1 + pred2
formula <- counts ~ 1

##FPBKpred(formula = formula, data = data, xcoordcol = xcoordcol,
##  ycoordcol = ycoordcol, covstruct = "Exponential", FPBK.col = NULL)
