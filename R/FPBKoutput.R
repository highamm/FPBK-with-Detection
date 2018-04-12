#' Create maps and summaries from FPBK results.
#'
#' The main purpose of this function is to take the results from FPBK and make
#' readable maps, the fitted variogram model, and brief summaries.
#'
#' @param pred.info is the output from \code{FPBK.pred} in this package.
#' @param conf.level is the desired confidence level for the prediction
#' @return maps of the kriged and observed counts, the fitted model variogram,
#' and any requested normal-based confidence intervals.
#' @export FPBKoutput



FPBKoutput <- function(pred.info, conf.level = 0.95) {

pred.total <- pred.info[[1]]
pred.total.var <- pred.info[[2]]
pred.vals <- pred.info[[3]]
covparmests <- pred.info[[4]]

confbounds <- matrix(round(as.numeric(pred.total) + c(1, -1) *
    stats::qnorm((1 - conf.level) / 2) *sqrt(as.numeric(pred.total.var))), nrow = 1)

labs <- c("Lower Bound", "Upper Bound")
colnames(confbounds) <- labs

## still deciding what kinds of plots should be output....
##plot.fun <- function() {
##plot.df <- as.data.frame(pred.vals)
##ggplot2::ggplot(data = plot.df, ggplot2::aes(x = V1, y = V2, colour = V3, shape = as.factor(V4))) +
##  ggplot2::geom_point()

##x <- seq(0, 1, length.out = 3)
##y <- seq(0, 1, length.out = 5)
##coords <- expand.grid(x, y)
##vals <- rnorm(n = nrow(coords), mean = 5, sd = 3)
##test.df <- as.data.frame(cbind(coords, vals))
##test.df$Var3
##ggplot(data = test.df, aes(x = Var1, y = Var2)) + geom_tile(aes(fill = vals)) +
##  xlim(c(0, 2)) + ylim(c(0, 2))
##}
}

##pred.info <- FPBKpred(formula = formula, data = data, xcoordcol = xcoordcol,
##  ycoordcol = ycoordcol, covstruct = "Exponential", FPBK.col = NULL)
conf.level <- 0.95
