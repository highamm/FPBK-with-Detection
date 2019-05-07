#' Extract Model Residuals
#'
#' @param object a slmfit object
#' @param ... additional arguments
#' @return a vector of residuals. Without detection, the residuals are  each observed count minus the estimated mean. With detection, the residuals are each observed count minus (estimated mean times the estimated detection).
#' @export

residuals.slmfit <- function(object, ...)
{
  resid.vec <- object$resids
  return(resid.vec)
}