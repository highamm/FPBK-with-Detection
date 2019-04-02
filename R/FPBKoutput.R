#' Create maps and summaries from FPBK results.
#'
#' The main purpose of this function is to take the results from FPBK and make
#' readable maps, a fitted variogram plot, and normal-based prediction intervals. The main input for this function is the output from the \code{FPBK.pred} function.
#'
#' @param pred_info is the output from \code{FPBK.pred} in this package.
#' @param conf_level is the desired confidence level for the prediction. If \code{conf_level} is a vector, then confidence intervals for 
#' each element of the vector will be produced.
#' @param get_krigmap is an indicator for whether or not a grid of
#' the kriged responses is returned
#' @param get_sampdetails is an indicator for whether or not a summary
#' of the sampled counts should be output. This summary includes 
#' the total number of animals or plants sited, the total area 
#' surveyed, the number of sampled sites, the total number of sites,
#' etc.
#' @param get_variogram is an indicator for whether or not
#' a variogram of the residuals should be returned
#' @param nbreaks is the number of breakpoints used in the spatial graphic
#' @param breakMethod is either \code{'quantile'} or \code{'even'} and determines how the break points are constructed.
#' @param pointsize is the size of the points on the spatial graphic.
#' @return \itemize{
#'   \item prediction for the total with prediction intervals
#'   \item a map of the kriged counts (optional)
#'   \item a summary of the sample data (optional)
#'   \item an empirical variogram from \code{gstat} with the fitted variogram model with details of the empirical variogram and spatial parameter estimates (optional)
#' }
#' @import ggplot2
#' @export FPBKoutput

## next steps: write up delta method
## consider scaling the nugget variance for the stratified sites
## add option to input shapefile
## make sure maps show counts, not densities

FPBKoutput <- function(pred_info, conf_level = c(0.80, 
  0.90, 0.95),
  get_krigmap = FALSE, get_sampdetails = FALSE,
  get_variogram = FALSE,
  nbreaks = 4,
  breakMethod = 'quantile', 
  pointsize = 2) {

pred.total <- pred_info$FPBK_Prediction
pred.total.var <- pred_info$PredVar
pred.vals <- data.frame(pred_info$Pred_df)
covparmests <- pred_info$SpatialParms
formula <- pred_info$formula
respname <- base::all.vars(formula)[1]
sampind <- pred.vals[ ,paste(base::all.vars(formula)[1], "_sampind",
  sep = "")]
preds <- pred.vals[ ,paste(base::all.vars(formula)[1], "_pred",
  sep = "")]
responsevar <- pred.vals[ ,base::all.vars(formula)[1]]
CorModel <- pred_info$CorModel
muhat <- pred.vals[ ,paste(base::all.vars(formula)[1], "_muhat",
  sep = "")]
piest <- pred.vals[ ,paste(base::all.vars(formula)[1], "_piest",
  sep = "")]
confbounds <- matrix(NA, nrow = length(conf_level), ncol = 3)
areavar <- pred.vals[ ,paste(base::all.vars(formula)[1], "_areas",
  sep = "")]

simptab <- t(matrix(c(pred.total, sqrt(pred.total.var))))
colnames(simptab) <- c("Prediction", "SE(Prediction)")
print(simptab)

for (k in 1:length(conf_level)){
confbounds[k, ] <- matrix(c(round(as.numeric(pred.total) + c(1, -1) *
    stats::qnorm((1 - conf_level[k]) / 2) *
    sqrt(as.numeric(pred.total.var))), 
  as.vector(round(stats::qnorm((1 - conf_level[k]) / 2) * -1 * 
      sqrt(as.numeric(pred.total.var)) / pred.total, 2))), nrow = 1)
}

labs <- c("Lower Bound", "Upper Bound", "Proportion of Mean")

rowlabs <- rep(NA, length(conf_level))
for (j in 1:length(conf_level)) {
  rowlabs[j] <- paste(conf_level[j] * 100, "%")
}

colnames(confbounds) <- labs
rownames(confbounds) <- rowlabs
print(confbounds)

basicpred <- t(matrix(round(c(pred.total, sqrt(pred.total.var)))))
colnames(basicpred) <- c("Predicted Total", "SE(Total)")

if (get_sampdetails == TRUE) {
  nsitessampled <- sum(sampind)
  nsitestotal <- nrow(pred.vals)
  animalscounted <- sum(responsevar[sampind == 1])
  totalareacol <- sum(areavar)
  totalarea <- sum(areavar)
  areasampled <- sum(areavar[sampind == 1])
  outptmat <- t(matrix(c(nsitessampled, nsitestotal, animalscounted,
    totalarea, areasampled)))
  colnames(outptmat) <- c("Numb. Sites Sampled", "Total Numb. Sites",
    "Animals Counted", "Total Area", "Area Sampled")
  print(outptmat)
  
} else{
    outptmat <- NULL
  }

if (get_variogram == TRUE) {
  sampled_df <- data.frame(subset(pred.vals, sampind == 1))
  muhatsamp <- muhat[sampind == 1]
  predsamp <- preds[sampind == 1] / areavar[sampind == 1]
  piestsamp <- piest[sampind == 1]
  responsevarsamp <- responsevar[sampind == 1]
  sampled_df$resids <- as.vector(responsevarsamp / (piestsamp *
      areavar[sampind == 1]) -
      muhatsamp / piestsamp)
  ##sampled_df$resids <- as.vector(predsamp -
  ##    muhatsamp / piestsamp)
  ## not sure what to do about residuals here since we do not
  ## actually observe any totals
  
  ## code for empirical variogram
  g_obj <- gstat::gstat(formula = resids ~ 1,
    locations = ~ xcoordsUTM_ + ycoordsUTM_,
    data = sampled_df)
  vario_out <- gstat::variogram(g_obj, cressie = FALSE)##,
    ##cutoff = 150)
  maxy <- max(vario_out$gamma)
  
  vartab <- cbind(vario_out$dist, vario_out$gamma, 
    vario_out$np)
  colnames(vartab) <- c("Distance", "Gamma", "Number of Pairs")
  covparmmat <- t(matrix(covparmests))
  colnames(covparmmat) <- c("Nugget", "Partial Sill", "Range")
  ## code for fitted variogram
  
  maxdist <- max(vario_out$dist)
  x.dist.plot <- seq(0, maxdist, length.out = 300)
  nugget <- covparmests[1]
  
  if (CorModel == "Exponential") {
  v.modfit <- nugget + covparmests[2] -
    covparmests[2] * corModelExponential(x.dist.plot, covparmests[3])
  } else if (CorModel == "Gaussian") {
  v.modfit <- nugget + covparmests[2] -
      covparmests[2] * corModelGaussian(x.dist.plot, covparmests[3])
  } else if (CorModel == "Spherical") {
  v.modfit <- nugget + covparmests[2] - 
    covparmests[2] * corModelSpherical(x.dist.plot, covparmests[3])
  }
  tab2 <- cbind(x.dist.plot, v.modfit)
  df.plot <- as.data.frame(tab2)

  plot_out <- ggplot(data = vario_out,
    aes_(x = ~dist, y = ~gamma)) +
    geom_point(aes_(size = ~gstat::variogram(g_obj)$np)) +
    ylim(0, maxy * (15 / 14)) +
    geom_line(data = df.plot, aes_(x = ~x.dist.plot, y = ~v.modfit)) +
    xlab("Distance (UTM)") +
    ylab("Semi-Variogram") +
    ggtitle(paste("Empirical Variogram with Fitted",
      CorModel, "Model")) +
    scale_size_continuous("Number of Pairs")
  
  print(plot_out)
} else if (get_variogram == FALSE) {
  plot_out <- NULL
  covparmmat <- NULL
  vartab <- NULL
}

if (get_krigmap == TRUE) {
  ## need to make some decisions about how complex the grid/
  ## map should be
  alldata <- data.frame(pred.vals)
  shapevals <- c(1, 16)
  # (p1 <- ggplot(data = alldata, aes_(x = ~xcoords, y = ~ycoords,
  #   colour = ~preds, shape = ~as.factor(sampind))) + 
  #   geom_point(size = 4) +
  #   scale_colour_gradient2(low = "blue", mid = "yellow", high = "red",
  #     midpoint = median(alldata$preds)) +
  #     scale_shape_manual(values = shapevals))

  ## make rectangles based on minimum distance
  ## this will only work if data form a grid

  #minxdist <- min(dist(alldata$xcoords)[dist(alldata$xcoords) != 0])
  #minydist <- min(dist(alldata$ycoords)[dist(alldata$ycoords) != 0])
 nameresp <- as.vector((noquote(paste(base::all.vars(formula)[1],
   "_pred",
    sep = ""))))
 alldata$sampindfact_ <- factor(sampind)
 
 # name of column with predictions
 pcolname = paste(base::all.vars(formula)[1], "_pred",
   sep = "")

 # create breaks according to user-specified breakMethod
 
 if(breakMethod == 'quantile') {
   probs = (1:nbreaks)/(nbreaks + 1)
   brks = c(min(pred.vals[ ,pcolname]) - 1e-10,
     quantile(pred.vals[ ,pcolname], probs = probs),
     max(pred.vals[ ,pcolname]) + 1e-10)
 }
 if(breakMethod == 'even') {
   rang = max(pred.vals[ ,pcolname]) + 1e-10 -
     min(pred.vals[ ,pcolname]) - 1e-10
   brks = c(min(pred.vals[ ,pcolname]) - 1e-10,
     min(pred.vals[ ,pcolname]) - 1e-10 +
       rang*(1:nbreaks)/(nbreaks + 1),
     max(pred.vals[ ,pcolname]) + 1e-10)
 }
 # cut predictions at breakpoints to create a vector of factors and labels
 cuts = cut(pred.vals[,pcolname], breaks = brks)
 # create a color palette
## palette = viridisLite::viridis(length(levels(cuts)))

 
 p3 <- ggplot2::ggplot(data = alldata, aes_(x = ~xcoordsUTM_,
   y = ~ycoordsUTM_, shape = ~sampindfact_)) +  ##)) + 
   geom_point(aes(colour = cuts), size = pointsize,
     stroke = 3) +
   scale_fill_viridis_d() +
   scale_colour_viridis_d(name = "Counts") + 
   theme_bw() +
   scale_shape_manual("Samp Indicator", 
     labels = c("Unsampled", "Sampled"), values = shapevals) +
   theme(panel.background = element_blank(),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank()) +
   xlab("") + ylab("")
 
 
 
 
 
 # ## change this for binning like Jay did in the other package!
 #  p3 <- ggplot2::ggplot(data = alldata, aes_(x = ~xcoordsUTM_, y = ~ycoordsUTM_, shape = ~sampindfact_)) +  ##)) + 
 #    geom_point(aes_string(colour = nameresp)) +
 #     ## geom_rect(aes_(xmin = ~ (xcoords - minxdist / 2),
 #      ##  xmax = ~(xcoords + minxdist / 2),
 #      ##  ymin = ~(ycoords - minydist / 2),
 #      ##  ymax = ~(ycoords + minydist / 2),
 #      ##  fill = ~preds)) +
 #      scale_fill_viridis_c() +
 #      scale_colour_viridis_c() +
 #      scale_shape_manual(values = shapevals)
  
  print(p3)
  
} else {
  p3 <- NULL
}

tabsandfigs <- list(simptab, confbounds, outptmat, plot_out,
  covparmmat, vartab, p3)
names(tabsandfigs) <- c("basic", "conf", "suminfo",
  "varplot", "covparms",
  "varplottab", "krigmap")
return(tabsandfigs)
}

##pred_info <- FPBKpred(formula = formula, data = data, xcoordcol = xcoordcol,
##   ycoordcol = ycoordcol, CorModel = "Gaussian",
##  coordtype = "UTM", areacol = "areavar", FPBKcol = NULL)
# pred_info$Pred_df[ ,5]
# 
# FPBKoutput(pred_info = pred_info, get_variogram = TRUE)
# conf_level <- 0.95
