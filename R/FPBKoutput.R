#' Create maps and summaries from FPBK results.
#'
#' The main purpose of this function is to take the results from FPBK and make
#' readable maps, the fitted variogram model, and brief summaries.
#'
#' @param pred_info is the output from \code{FPBK.pred} in this package.
#' @param conf_level is the desired confidence level for the prediction
#' @param get_krigmap is an indicator for whether or not a grid of
#' the kriged responses is returned
#' @param get_variogram is an indicator for whether or not
#' a variogram of the residuals should be returned
#' @param CorModel Covariance model used, which is required to obtain 
#' the appropriate fitted model semi-variogram.
#' @return maps of the kriged and observed counts, the fitted model variogram,
#' and any requested normal-based confidence intervals.
#' @import ggplot2
#' @export FPBKoutput

## next steps: write up delta method
## consider scaling the nugget variance for the stratified sites
## add option to input shapefile
## make sure maps show counts, not densities

FPBKoutput <- function(pred_info, conf_level = 0.95,
  get_krigmap = FALSE,
  get_variogram = FALSE, CorModel = "Exponential") {

pred.total <- pred_info[[1]]
pred.total.var <- pred_info[[2]]
pred.vals <- data.frame(pred_info[[3]])
covparmests <- pred_info[[4]]

confbounds <- matrix(round(as.numeric(pred.total) + c(1, -1) *
    stats::qnorm((1 - conf_level) / 2) *
    sqrt(as.numeric(pred.total.var))), nrow = 1)

labs <- c("Lower Bound", "Upper Bound")
colnames(confbounds) <- labs

if (get_variogram == TRUE) {
  sampled_df <- data.frame(subset(pred.vals, pred.vals[ ,4] == 1))
  sampled_df$resids <- as.vector(sampled_df$preds - sampled_df$muhat)

  ## code for empirical variogram
  g_obj <- gstat::gstat(formula = resids ~ 1,
    locations = ~ xcoords + ycoords,
    data = sampled_df)
  vario_out <- gstat::variogram(g_obj)
  maxy <- max(vario_out$gamma)
  
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

  plot_out <- ggplot2::ggplot(data = vario_out,
    aes(x = dist, y = gamma)) +
    geom_point(aes(size = gstat::variogram(g_obj)$np)) +
    ylim(0, maxy * (15 / 14)) +
    geom_line(data = df.plot, aes(x = x.dist.plot, y = v.modfit)) +
    xlab("Distance (UTM)") +
    ylab("Semi-Variogram") +
    ggtitle(paste("Empirical Variogram with Fitted",
      CorModel, "Model")) +
    scale_size_continuous("Number of Pairs")
  
  print(plot_out)
}

if (get_krigmap == TRUE) {
  ## need to make some decisions about how complex the grid/
  ## map should be
  alldata <- data.frame(pred.vals)
  shapevals <- c(16, 15)
  (p1 <- ggplot(data = alldata, aes_(x = ~xcoords, y = ~ycoords,
    colour = ~preds, shape = ~as.factor(sampind))) + 
    geom_point(size = 4) +
    scale_colour_gradient2(low = "blue", mid = "yellow", high = "red",
      midpoint = median(alldata$preds)) +
      scale_shape_manual(values = shapevals))

  ## make rectangles based on minimum distance
  ## this will only work if data form a grid

  minxdist <- min(dist(alldata$xcoords)[dist(alldata$xcoords) != 0])
  minydist <- min(dist(alldata$ycoords)[dist(alldata$ycoords) != 0])

  (p3 <- ggplot(data = alldata, aes_(x = ~xcoords, y = ~ycoords,
    colour = ~preds, shape = ~as.factor(sampind))) + 
      geom_rect(aes_(xmin = ~ (xcoords - minxdist / 2),
        xmax = ~(xcoords + minxdist / 2),
        ymin = ~(ycoords - minydist / 2),
        ymax = ~(ycoords + minydist / 2),
        fill = ~preds)) +
      scale_fill_viridis_c() +
      scale_colour_viridis_c() +
      scale_shape_manual(values = shapevals))
  
  print(p3)
  
}


}

##pred_info <- FPBKpred(formula = formula, data = data, xcoordcol = xcoordcol,
##   ycoordcol = ycoordcol, CorModel = "Gaussian",
##  coordtype = "UTM", areacol = "areavar", FPBKcol = NULL)
# pred_info$Pred_df[ ,5]
# 
# FPBKoutput(pred_info = pred_info, get_variogram = TRUE)
# conf_level <- 0.95
