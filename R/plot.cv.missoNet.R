#' Plot the cross-validation errors produced by cv.missoNet
#' 
#' S3 method for plotting the cross-validation error surface from a fitted \code{'cv.missoNet'} object.
#' 
#' @param x A fitted \code{'cv.missoNet'} object.
#' @param type A character string for the type of plot, can be either "\code{cv.heatmap}" (default) or "\code{cv.scatter}".
#' @param detailed.axes Logical: whether the detailed axes should be plotted. The default is \code{'TRUE'}.
#' @param plt.surf Logical: whether to draw the error surface. The default is \code{'TRUE'}. This is only needed when \code{'type'} = "\code{cv.scatter}".
#' @param ... Other graphical arguments used by \sQuote{ComplexHeatmap::Heatmap} (\code{'type'} = "\code{cv.heatmap}") or \sQuote{scatterplot3d::scatterplot3d} (\code{'type'} = "\code{cv.scatter}").
#' 
#' @return The plot object.
#' @method plot cv.missoNet
#' @export
#' 
#' @author Yixiao Zeng \email{yixiao.zeng@@mail.mcgill.ca}, Celia M.T. Greenwood and Archer Yi Yang.
#' 
#' @examples
#' ## Simulate a dataset.
#' set.seed(123)  # reproducibility
#' sim.dat <- generateData(n = 200, p = 10, q = 10, rho = 0.1, missing.type = "MCAR")
#' 
#' \donttest{
#' ## Perform a five-fold cross-validation on the simulated dataset.
#' cvfit <- cv.missoNet(X = sim.dat$X, Y = sim.dat$Z, kfold = 5,
#'                      fit.1se = TRUE, permute = TRUE, with.seed = 486)
#' 
#' 
#' ## Plot the (standardized) mean cross-validated errors in a heatmap.
#' plot(cvfit, type = "cv.heatmap")
#' 
#' ## Plot the (standardized) mean cross-validated errors in a 3D scatterplot.
#' plot(cvfit, type = "cv.scatter", plt.surf = TRUE)
#' }

plot.cv.missoNet <- function(x, type = c("cv.heatmap", "cv.scatter"), detailed.axes = TRUE, plt.surf = TRUE, ...) {
  type <- match.arg(type)
  switch(type, cv.heatmap = plot.heatmap(x, detailed.axes, ...), cv.scatter = plot.scatter(x, detailed.axes, plt.surf, ...))
}
