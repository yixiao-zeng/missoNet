#' Plot the cross-validation surface produced by cv.missoNet
#' @param x Fitted cv.missoNet object.
#' @param type Type of plot, can be either "cv.heatmap" or "cv.scatter". The default is "cv.heatmap".
#' @param ... Other graphical parameters.
#' @method plot cv.missoNet
#' @export
#' @examples 

plot.cv.missoNet <- function(x, type = c("cv.heatmap", "cv.scatter"), ...) {
  type <- match.arg(type)
  switch(type, cv.heatmap = plot.heatmap(x), cv.scatter = plot.scatter(x))
}
