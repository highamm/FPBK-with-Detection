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
#' @param detinfo is a vector consisting of the mean detection
#' probability and its standard error. By default, it is set to
#' \code{c(1, 0)}, indicating perfect detection (1) with no
#' variance (0). Using the default assumes perfect detection.
#' @param areacol is the name of the column with the areas of the sites. By default, we assume that all sites have equal area, in which
#' case a vector of 1's is used as the areas.
#' @param FPBKcol is the name of the column with the prediction 
#' weights if the analyst would like to predict something other
#' than the population total, like the total for a subset of the
#' entire region of interest.
#' @param stratcol is the column in the data set that contains the stratification variable. 
#' @return a report with information about the predicted total across all sites as well as variogram information for each stratum in the column \code{stratcol}.
#' @import stats
#' @export multistrat


multistrat <- function(formula, data, xcoordcol, ycoordcol,
  CorModel = "Exponential",
  coordtype = "LatLon", estmethod = "REML",
  covestimates = c(NA, NA, NA),
  detectionobj = NULL,
  detinfo = c(1, 0), 
  areacol = NULL,
  FPBKcol = NULL,
  stratcol = NULL) {
  
  stratvar <- factor(data[ ,stratcol])

  slmfitouts <- vector("list",  length(levels(stratvar)))
  predictouts <- vector("list",  length(levels(stratvar)))
  dfout <- vector("list", length(levels(stratvar)))
  prediction <- vector("list", length(levels(stratvar)))
  predictionvar <- vector("list", length(levels(stratvar)))
  tabsout <- vector("list", length(levels(stratvar)))
  muhat <- vector("list", length(levels(stratvar)))
  corrpimat <- vector("list", length(levels(stratvar)))
  tlambda <- vector("list", length(levels(stratvar)))
  maxy <- rep(NA, length(levels(stratvar)))
  
  for (k in 1:length(levels(stratvar))) {
  slmfitouts[[k]] <- slmfit(formula = formula,
    data = subset(data, stratvar == levels(stratvar)[k]),
    xcoordcol = xcoordcol, ycoordcol = ycoordcol,
    CorModel = CorModel, 
    coordtype = coordtype,
     estmethod = estmethod,
    covestimates = covestimates, detectionobj = detectionobj,
    areacol = areacol)
  
  predictouts[[k]] <- predict(object = slmfitouts[[k]],
    FPBKcol = FPBKcol, detinfo = detinfo)
  
  if (is.null(detectionobj) == FALSE) {
    sampind <- predictouts[[k]]$Pred_df[ ,paste(base::all.vars(formula)[1], "_sampind", sep = "")] 
    
    corrpimat[[k]] <- slmfitouts[[k]]$piboot[sampind, ]
    
    muhat[[k]] <- predictouts[[k]]$Pred_df[sampind , paste(base::all.vars(formula)[1], "_muhat", sep = "")]
    
    tlambda[[k]] <- predictouts[[k]]$tlambda
    
  } else {
    
  }
  

  ## essentially have to store the report for each k,
  ## then append these reports to the final report for all
  ## of the strata that is generated outside of the loop
  tabsout[[k]] <- FPBKoutput(predictouts[[k]], conf_level = c(0.80, 
    0.90, 0.95),
    get_krigmap = FALSE, get_sampdetails = TRUE,
    get_variogram = TRUE,
   ## nbreaks = 4,
  ##  breakMethod = 'quantile', 
    pointsize = 2)
  stratname <- levels(stratvar)[k]
  rownames(tabsout[[k]]$basic) <- stratname
  rownames(tabsout[[k]]$conf) <- paste(stratname,
    rownames(tabsout[[k]]$conf))
  rownames(tabsout[[k]]$suminfo) <- stratname
  rownames(tabsout[[k]]$covparms) <- stratname
  rownames(tabsout[[k]]$varplottab) <- rep(stratname, 
    nrow(tabsout[[k]]$varplottab))
  tabsout[[k]]$varplot <- tabsout[[k]]$varplot + ggplot2::ylab(paste("Semi-Variogram", stratname))

  maxy[k] <- (7 / 6) * max(tabsout[[k]]$varplottab[ ,2])

  prediction[[k]] <- predictouts[[k]]$FPBK_Prediction
  predictionvar[[k]] <- predictouts[[k]]$PredVar
  dfout[[k]] <- predictouts[[k]]$Pred_df
  ## could have a page for the low stratum report, a page for the
  ## high stratum report, and then a combined report....it will take 
  ## a bit of work combining the results but not too much I think
  }
  
  ##stratvartest <- factor(c("A", "B", "C"))
  
  ## next: incorporate this into the FPBKoutput only for
  ## this multistrat function part
  
  if (is.null(detectionobj) == FALSE) {
    
      if(detectionobj$varmethod == "Bootstrap") {
    extravar <- matrix(0, nrow = nlevels(stratvar),
      ncol = nlevels(stratvar))
    
    for (j in 1:nlevels(stratvar)) {
      for (jj in 1:nlevels(stratvar)) {
        if (j != jj) {
          extravar[j, jj] <- matrixcalc::hadamard.prod(tlambda[[j]], t(matrix(muhat[[j]]))) %*%
            cov(t(corrpimat[[j]]), t(corrpimat[[jj]])) %*%
            t(matrixcalc::hadamard.prod(tlambda[[jj]], t(matrix(muhat[[jj]]))))
        } else if (j == jj) {
          extravar[j, jj] <- 0
        }
      }
    }
  } else if(detectionobj$varmethod == "Delta") {
    extravar <- 0
    warning("The confidence interval from using the delta method to obtain the variance of the detection probabilities will be slightly too small due to not accounting for correlation across strata due to using the same detection data.")
  } 
    } else {
    extravar <- 0
  }
  
  
  maxall <- max(maxy)
  
  allFPBKobj <- vector("list", length(predictouts[[k]]))
  allFPBKobj[[1]] <- matrix(do.call("sum", prediction))
  ## need to update the variance to account for covariance
  ## from using the same detection for each stratum
  allFPBKobj[[2]] <- matrix(do.call("sum", predictionvar) +
      sum(extravar)) 
  allFPBKobj[[3]] <- do.call("rbind", dfout)
  allFPBKobj[[4]] <- c(NA, NA, NA)
  allFPBKobj[[5]] <- formula
  allFPBKobj[[6]] <- "DontNeed"

  names(allFPBKobj) <- c("FPBK_Prediction", "PredVar",
    "Pred_df", "SpatialParms", "formula", "CorModel")
  
  tabsall <- FPBKoutput(allFPBKobj, conf_level = c(0.80, 
    0.90, 0.95),
    get_krigmap = TRUE, get_sampdetails = TRUE,
    get_variogram = FALSE,
  ##  nbreaks = 4,
  ##  breakMethod = 'quantile', 
    pointsize = 7)

  rownames(tabsall$basic) <- "Total"
  rownames(tabsall$conf) <- rep("Total", nrow(tabsout[[1]]$conf))
  rownames(tabsall$suminfo) <- "Total"
  tabsall$krigmap <- tabsall$krigmap + ggplot2::ggtitle("Map of Predictions for All Sites")


  basicinfotab <- vector("list", nlevels(stratvar))
  basicinfotab[[1]] <- round(tabsout[[1]][[1]], 2)
  conftab <- vector("list", nlevels(stratvar))
  conftab[[1]] <- tabsout[[1]][[2]]
  suminform <- vector("list", nlevels(stratvar))
  suminform[[1]] <- round(tabsout[[1]][[3]], 0)
  covparmtab <- vector("list", nlevels(stratvar))
  covparmtab[[1]] <- round(tabsout[[1]]$covparms, 3)
  vartab <- vector("list", nlevels(stratvar))
  vartab[[1]] <- round(tabsout[[1]]$varplottab, 2)
  
    for (k in 2:length(levels(stratvar))) {
      basicinfotab[[k]] <- rbind(basicinfotab[[k - 1]],
        round(tabsout[[k]][[1]], 2))
      conftab[[k]] <- rbind(conftab[[k - 1]], 
        rep(".", ncol(tabsout[[1]][[2]])),
        tabsout[[k]][[2]])

      suminform[[k]] <- rbind(suminform[[k - 1]],
        round(tabsout[[k]][[3]], 0))
      
      covparmtab[[k]] <- rbind(covparmtab[[k - 1]],
        round(tabsout[[k]]$covparms, 3))
      
     vartab[[k]] <- rbind(vartab[[k - 1]], rep(".", ncol(tabsout[[1]]$varplottab)), rep(".", ncol(tabsout[[1]]$varplottab)),
        round(tabsout[[k]]$varplottab, 2))
     
     dfout[[k]] <- rbind(dfout[[k - 1]], dfout[[k]])
     
      ## keep variogram tables separate
    }
  

  varplots <- vector("list", nlevels(stratvar))
  for (k in 1:nlevels(stratvar)) {
 varplots[[k]] <- suppressMessages(tabsout[[k]]$varplot + ggplot2::ylim(c(0, maxall)))
  }
  
  
  basic <- rbind(basicinfotab[[k]], round(tabsall$basic, 2))
  conf <- rbind(conftab[[k]],
    rep(".", ncol(conftab[[k]])), tabsall$conf)
  suminfo <- rbind(suminform[[k]], round(tabsall$suminfo, 0))
  varplot <- varplots
    ##varplot <- gridExtra::grid.arrange(grobs = varplots, nrow = 1)
  covparms <- covparmtab[[k]]
  varplottab <- vartab[[k]]
  krigmap <- tabsall$krigmap
  pred.vals <- dfout[[k]]
  
  output_info <- list(basic, conf, suminfo, varplot, 
    covparms, varplottab, krigmap, pred.vals)
  names(output_info) <- c("basic", "conf",
    "suminfo", "varplot", "covparms", "varplottab",
    "krigmap", "predvals")
  return(output_info)
  #get_reportdoc(output_info = output_info)
  
  ## export all figures into their own jpeg files
  ## keep as html...can change to PDF using print
  
 ## ggplot2::scale_shape_manual("Samp Indicator", 
  ##    labels = c("Unsampled", "Sampled", "x", "xx"),
  ##    values = shapevals)

}


