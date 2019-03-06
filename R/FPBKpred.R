#' Perform Finite Population Block Kriging
#'
#' Takes the output of \code{slmfit} and uses FPBK to predict the counts on the unsampled sites.
#' The column with the counts should have numeric values for the observed counts
#' on the sampled sites and `NA` for any site that was not sampled.
#'
#'
#' @param object is an object generated from \code{slmfit}
#' @param FPBKcol is the name of the column that contains the weights for prediction. The default setting predicts the population total.
#' @param detinfo is a vector consisting of the mean detection
#' probability and its standard error. By default, it is set to
#' \code{c(1, 0)}, indicating perfect detection (1) with no
#' variance (0). Using the default assumes perfect detection.
#' @param ... Additional arguments
#' @return a list with \itemize{
#'   \item the estimated population total
#'   \item the estimated prediction variance
#'   \item a data frame containing \enumerate{
#'        \item x-coordinates
#'        \item y-coordinates
#'        \item density predictions
#'        \item indicator variable for whether or not the each site was sampled
#'        \item estimated mean for each site
#'        }
#'    \item vector with estimated covariance parameters
#' }
#' @import stats
#' @export


predict.slmfit <- function(object, FPBKcol = NULL,
  detinfo = c(1, 0), ...) {
  
  ## if FPBKcol is left out, we are predicting the population total.
  ## Otherwise, FPBKcol is the name of the column in the data set
  ## with the weights for the sites that we are predicting (eg. a vector
  ## of 1's and 0's for predicting the total of the sites marked with 1's)
  
  
  formula <- object$FPBKpredobj$formula
  data <- object$FPBKpredobj$data
  xcoordsUTM <- object$FPBKpredobj$xcoordsUTM
  ycoordsUTM <- object$FPBKpredobj$ycoordsUTM
  covparmests <- object$SpatialParmEsts
  
  if (is.null(FPBKcol) == TRUE) {
    predwts <- rep(1, nrow(data))
  } else if (is.character(FPBKcol) == TRUE) {
    predwts <- data[ ,FPBKcol]
  } else{
    stop("FPBKcol must be a character specifying the name of the
      column of
      prediction weights in the data set")
  }
  
  
  
  
  
  
  fullmf <- stats::model.frame(formula, na.action =
      stats::na.pass, data = data)
  yvar <- stats::model.response(fullmf, "numeric")
  density <- yvar
  
  ind.sa <- !is.na(yvar)
  ind.un <- is.na(yvar)
  data.sa <- data[ind.sa, ]
  data.un <- data[ind.un, ]
  
  B <- predwts
  Bs <- B[ind.sa]
  Bu <- B[ind.un]
  
  m.un <- stats::model.frame(formula, data.un, na.action =
      stats::na.pass)
  Xu <- stats::model.matrix(formula, m.un)
  X <- stats::model.matrix(formula, fullmf)
  
  ## sampled response values and design matrix
  m.sa <- stats::model.frame(formula, data.sa, na.action =
      stats::na.omit)
  z.sa <- stats::model.response(m.sa)
  Xs <- stats::model.matrix(formula, m.sa)
  z.density <- z.sa
  
  if (object$DetectionInd == TRUE) {
  
  
  Sigma <- object$FPBKpredobj$covmat
  
  ## used in the Kriging formulas
  Sigma.us <- Sigma[ind.un, ind.sa]
  Sigma.su <- t(Sigma.us)
  Sigma.ss <- Sigma[ind.sa, ind.sa]
  Sigma.uu <- Sigma[ind.un, ind.un]
  
  ## give warning if covariance matrix cannot be inverted
  # if(abs(det(Sigma.ss)) <= 1e-21) {
  #   warning("Covariance matrix is compulationally singular and
  #     cannot be inverted")
  # }
  
  Sigma.ssi <- object$FPBKpredobj$covmatsampi
  
  ## the generalized least squares regression coefficient estimates
  betahat <- object$CoefficientEsts
  
  ## estimator for the mean vector
  muhats <- (Xs %*% betahat) 
  muhatu <- (Xu %*% betahat) 
  
  muhat <- rep(NA, nrow(data))
  muhat[ind.sa == TRUE] <- muhats / detinfo[1]
  muhat[ind.sa == FALSE] <- muhatu / detinfo[1]
  
  
  ## matrices used in the kriging equations
  ## notation follow Ver Hoef (2008)
  Cmat <- Sigma.ss %*% as.matrix(Bs) +
    Sigma.su %*% as.matrix(Bu)
  Dmat <- t(X) %*% matrix(B) - t(Xs) %*% Sigma.ssi %*% Cmat
  Vmat <- solve(t(Xs) %*% Sigma.ssi %*% Xs)
  
  ## the predicted values for the sites that were not sampled
  zhatu <- (Sigma.us %*% Sigma.ssi %*% (z.density -
      muhats) + muhatu) / detinfo[1]
  
  ## creating a column in the outgoing data set for predicted counts as
  ## well as a column indicating whether or not the observation was sampled
  ## or predicted
  preddensity <- density / detinfo[1]
  preddensity[is.na(preddensity) == TRUE] <- zhatu
  
  sampind <- rep(1, length(yvar))
  sampind[is.na(yvar) == TRUE] <- 0
  
  ## adding the site-by-site predictions
  
  W <- (t(Xu) - t(Xs) %*% Sigma.ssi %*% Sigma.su) / detinfo[1]
  
  ## NEED TO CHANGE THIS FOR THE DETECTION CASE
  sitecov <- Sigma.uu - Sigma.us %*% Sigma.ssi %*% Sigma.su +
    t(W) %*% Vmat %*% W
  sitevarnodet <- diag(sitecov)
  
  
  densvar <- rep(NA, nrow(data))
  densvar[sampind == 1] <- 0
  densvar[sampind == 0] <- sitevarnodet
  
  varpibar <- detinfo[2] ^ 2
  varinvdelta <- varpibar * (1 / detinfo[1]) ^ 4
  
  sitevar <- (preddensity ^ 2) * varinvdelta +
    ((1 / detinfo[1]) ^ 2) * densvar +
    densvar * varinvdelta
  
  ## the FPBK predictor: already adjusted for detection
  FPBKpredictor <- (t(B) %*% preddensity)
  
  ## the prediction variance for the FPBK predictor
  ## if detectionest is left as the default, we assume
  ## perfect detection with no variability in the detection estimate
  
  pred.var <- (t(as.matrix(B)) %*% Sigma %*%
      as.matrix(B) -
      t(Cmat) %*% Sigma.ssi %*% Cmat +
      t(Dmat) %*% Vmat %*% Dmat)
  


  
  ## now adjusted for mean detection
  
  pred.var.obs <- (FPBKpredictor ^ 2) * varinvdelta +
    ((1 / detinfo[1]) ^ 2) * pred.var +
    pred.var * varinvdelta
  
  ## returns a list with 3 components:
  ## 1.) the kriging predictor and prediction variance
  ## 2.) a matrix with x and y coordinates, kriged predctions, and
  ## indicators for whether sites were sampled or not
  ## 3.) a vector of the estimated spatial parameters
  
  } else if (object$DetectionInd == FALSE) {
    
    R <- object$FPBKpredobj$R
    Xnstar <- object$FPBKpredobj$Xnstar
    Cssi <- object$FPBKpredobj$Cssi
    C <- object$FPBKpredobj$C
    D <- object$FPBKpredobj$covmat
    piest <- object$FPBKpredobj$piest
    
    ## breaking up calculation into a few different parts
    p10 <- t(B) %*% t(R) %*% Cssi %*% z.density

    p11 <- t(B) %*% X - t(B) %*% t(R) %*% Cssi %*% Xnstar
    
    p12 <- solve(t(Xnstar) %*% Cssi %*% Xnstar) %*%
      t(Xnstar) %*% Cssi %*% z.density
    
    FPBKpredictor <- p10 + p11 %*% p12
    tlambda <- t(B) %*% t(R) %*% Cssi + (t(B) %*% X
      - t(B) %*% t(R) %*% Cssi %*% Xnstar) %*%
      solve(t(Xnstar) %*% Cssi %*% Xnstar) %*% t(Xnstar) %*% Cssi
    
    pred.check <- tlambda %*% z.density

    p13 <- tlambda %*% C %*% t(tlambda)
    p14 <- t(tlambda %*% R %*% B)
    p15 <- tlambda %*% R %*% B
    p16 <- t(B) %*% D %*% B
    pred.var.obs <- p13 - p14 - p15 + p16
   ## se <- sqrt(varhat)
    
    ## stuff to be appended to the data frame
    ## estimator for the mean vector
    
    betahat <- object$CoefficientEsts
    
    muhatsz <- Xs %*% betahat
    muhatuz <- Xu %*% betahat
    
    muhats <- muhatsz * piest[ind.sa]
    muhatu <- muhatuz * piest[ind.un]
    muhat <- rep(NA, nrow(data))
    muhat[ind.sa == TRUE] <- muhats
    muhat[ind.sa == FALSE] <- muhatu
    
    sampind <- rep(1, length(yvar))
    sampind[is.na(yvar) == TRUE] <- 0
    
    
    weights <- t(R) %*% Cssi + (X - t(R) %*% Cssi %*% Xnstar) %*%
      solve(t(Xnstar) %*% Cssi %*% Xnstar) %*% t(Xnstar) %*% Cssi
    
    varweights <- weights %*% C %*% t(weights) -
      t(R) %*% t(weights) - weights %*% R + D
    ##varweightsarea <- varweights * (as.matrix(area.all) %*% t(as.matrix(area.all)))
    ##varhat.ML.bin.area <- t(B) %*% varweights %*% B
    
    
    pred.persite <- (weights %*% z.density)
    varweightsdiag <- diag(varweights)
    

    preddensity <- pred.persite
    densvar <- varweightsdiag

  } 
  
  
  df_out <- data.frame(cbind(data, xcoordsUTM, ycoordsUTM,
    preddensity, densvar, sampind, muhat))
  
  # data <- data.frame(y = 1:10, x = 2:11)
  #
  # fullmf <- stats::model.frame(formula, na.action =
  #   stats::na.pass, data = data)
  
  colnames(df_out) <- c(colnames(data), "_xcoordsUTM_", "_ycoordsUTM_",
    paste(base::all.vars(formula)[1], "_pred",
      sep = ""),
    paste(base::all.vars(formula)[1], "_predvar",
      sep = ""),
    paste(base::all.vars(formula)[1], "_sampind",
      sep = ""),
    paste(base::all.vars(formula)[1], "_muhat",
      sep = ""))
  obj <- list(FPBKpredictor, pred.var.obs,
    df_out,
    as.vector(covparmests),
    formula = formula)
  
  names(obj) <- c("FPBK_Prediction", "PredVar",
    "Pred_df", "SpatialParms", "formula")
  
  class(obj) <- "sptotalPredOut"
  
  return(obj)
  
  }
