#' Plot the cross-validation errors produced by cv.missoNet
#' 
#' S3 method for plotting the cross-validation error surface from a fitted \code{'cv.missoNet'} object.
#' 
#' @param x A fitted \code{'cv.missoNet'} object.
#' @param type A character string for the type of plot, can be either "\code{cv.heatmap}" (default) or "\code{cv.scatter}".
#' @param detailed.axes Logical: whether the detailed axes should be plotted. The default is \code{'TRUE'}.
#' @param ... Other graphical arguments used by \sQuote{ComplexHeatmap::Heatmap}.
#' @method plot cv.missoNet
#' @export
#' 
#' @author Yixiao Zeng \email{yixiao.zeng@@mail.mcgill.ca}, Celia M.T. Greenwood and Archer Yi Yang.
#' 
#' @examples
#' ## Perform a five-fold cross-validation on a simulated dataset.
#' sim.dat <- generateData(n = 200, p = 10, q = 10, rho = 0.1, missing.type = "MCAR")
#' cvfit <- cv.missoNet(X = sim.dat$X, Y = sim.dat$Z, kfold = 5, fit.1se = TRUE)
#' 
#' 
#' ## Plot the (standardized) mean cross-validated errors in a heatmap.
#' plot(cvfit)
#' 
#' 
#' ## Plot the (standardized) mean cross-validated errors in a scatterplot.
#' plot(cvfit, type = "cv.scatter", detailed.axes = FALSE)

plot.cv.missoNet <- function(x, type = c("cv.heatmap", "cv.scatter"), detailed.axes = TRUE, ...) {
  type <- match.arg(type)
  switch(type, cv.heatmap = plot.heatmap(x, detailed.axes, ...), cv.scatter = plot.scatter(x, detailed.axes))
}
