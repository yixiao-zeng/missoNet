#' Plot the cross-validation errors produced by \code{cv.missoNet}
#' 
#' S3 method for plotting the mean cross-validation error surface from a fitted `\code{cv.missoNet}` object.
#' 
#' @param x Fitted `\code{cv.missoNet}` object.
#' @param type Type of plot, can be either `\code{"cv.heatmap"}` (default) or `\code{"cv.scatter"}`.
#' @param detailed.axes Logical: should the detailed axes be plotted? The default is `TRUE`.
#' @param ... Other graphical arguments used by `\code{ComplexHeatmap::Heatmap}`.
#' @method plot cv.missoNet
#' @export

plot.cv.missoNet <- function(x, type = c("cv.heatmap", "cv.scatter"), detailed.axes = TRUE, ...) {
  type <- match.arg(type)
  switch(type, cv.heatmap = plot.heatmap(x, detailed.axes, ...), cv.scatter = plot.scatter(x, detailed.axes))
}
