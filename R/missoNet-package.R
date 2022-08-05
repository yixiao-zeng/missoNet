#' Multi-task regression and conditional network estimation with missing data
#'
#' \tabular{ll}{ Package: \tab missoNet\cr Type: \tab Package\cr Version: \tab
#' 0.0.0.9000\cr Date: \tab 2022-06-25\cr License: \tab GPL-2\cr }
#' The package 'missoNet' provides three categories of important functions:
#' cv.missoNet/missoNet; plot/predict and generateData.
#'
#' @section missoNet functions:
#' \describe{
#'   \item{cv.missoNet}{fit missoNet with the optimal regularization parameters selected 
#'     by a k-fold cross-validation over a grid of \code{lambda} values.}
#'   \item{missoNet}{fit missoNet with a given pair of regularization parameters.}
#'   \item{plot}{S3 method for plotting the cross-validation errors from a "cv.missoNet" object.}
#'   \item{predict}{S3 method for making predictions from a "cv.missoNet" object.}
#'   \item{generateData}{generate simulation data by specifying a certain type of missing mechanism 
#'     and missing probability.}
#' }
#'
#' @docType package
#' @name missoNet-package
#' @useDynLib missoNet, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics par
#' @importFrom stats coef complete.cases quantile rbinom rnorm runif sd
#' @importFrom utils install.packages setTxtProgressBar txtProgressBar
NULL
#> NULL