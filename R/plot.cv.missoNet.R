#' Plot the cross-validation errors produced by \code{cv.missoNet}
#' @param x Fitted \code{cv.missoNet} object.
#' @param type Type of plot, can be either \code{"cv.heatmap"} or \code{"cv.scatter"}. Default is \code{type = "cv.heatmap"}.
#' @param detailed.axis Logical: should the detailed axes be plotted? Default is \code{plot.axis = TRUE}.
#' @param ... Other graphical arguments used by \code{ComplexHeatmap}.
#' @method plot cv.missoNet
#' @export
#' @examples 

plot.cv.missoNet <- function(x, type = c("cv.heatmap", "cv.scatter"), detailed.axis = TRUE, ...) {
  type <- match.arg(type)
  switch(type, cv.heatmap = plot.heatmap(x, detailed.axis, ...), cv.scatter = plot.scatter(x, detailed.axis))
}
