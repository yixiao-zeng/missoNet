#' Multi-task regression and conditional network estimation with missing values in the tasks
#'
#' \tabular{ll}{ Package: \tab missoNet\cr Type: \tab Package\cr Version: \tab
#' 1.0.0\cr Date: \tab 2022-10-01\cr License: \tab GPL-2\cr }
#'
#' @section missoNet functions:
#' \describe{
#'   \item{missoNet}{Fit a series of \sQuote{\code{\link{missoNet}}} models with user-supplied regularization parameter pairs for the lasso penalties, \{(\eqn{\lambda_B}, \eqn{\lambda_\Theta})\}.}
#'   \item{cv.missoNet}{Perform k-fold cross-validation for \sQuote{\code{\link{missoNet}}} over a grid of (auto-computed) regularization parameter pairs.}
#'   \item{plot}{S3 method for plotting the cross-validation errors from a fitted \code{'cv.missoNet'} object.}
#'   \item{predict}{S3 method for making predictions of response values from a fitted \code{'cv.missoNet'} object.}
#'   \item{generateData}{Quickly generate synthetic data for simulation studies.}
#' }
#'
#' @docType package
#' @name missoNet-package
#' @useDynLib missoNet, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics par
#' @importFrom grid gpar grid.rect
#' @importFrom stats coef complete.cases quantile rbinom rnorm runif sd
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom glasso glasso
NULL
#> NULL