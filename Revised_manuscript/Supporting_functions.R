# function for taking the legend from the ggplot

#' g_legend
#'
#' @param a.gplot A plot that you would like to take the legend from
#'
#' @return
#' @export
#'
#' @examples
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#' get.st
#'
#' @param fit result of the optim function
#'
#' @return matrix with the values and upper and lower bands for the
#' @export
#'
#' @examples
get.st <- function(fit){
  # function to get the confidence interval on g0 and g1 given the output of optim
  fisher_info<-solve(-fit$hessian)
  prop_sigma<-sqrt(diag(fisher_info))

  return(prop_sigma)
}


