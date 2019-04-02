#' Function used to generate an HTML report.
#' 
#' This function produces an HTML report with the results from 
#' the finite population block kriging from \code{FPBKoutput}.
#'
#' @param output_info is the output from \code{FPBKoutput} in this package.
#' @import ggplot2
#' @export get_reportdoc

get_reportdoc <- function(output_info) {
  
file <- system.file("ReportTest.Rmd", package = "FPBKPack2")

## need to think more carefully about where this report
## should go.
dout <- "~/Desktop/"
if (missing(dout)) {
  dout <- getwd()
}

rmarkdown::render(file, envir = list(varplot = output_info$varplot, 
  varinfo = output_info$varplottab, krigplot = output_info$krigmap, 
  predtable = output_info$basic, conftable = output_info$conf,
  sumtable = output_info$suminfo, covparmests = output_info$covparms),
  output_dir = dout, 
  output_file = paste('report.', Sys.Date(), 
    '.html', sep=''))
}
