#' missoNet: a package for multivariate regression and conditional network estimation with missing outputs.
#'
#' The missoNet package provides three categories of important functions:
#' cv.missoNet/missoNet, plot/predict and generateData.
#' 
#' @section missoNet functions:
#' cv.missoNet: fit missoNet with the optimal regularization parameters selected by a k-fold cross-validation.
#' missoNet: fit missoNet with a given pair of regularization parameters.
#' plot: S3 method for plotting the cross-validation errors from a "cv.missoNet" object.
#' predict: S3 method for making predictions from a "cv.missoNet" object.
#' generateData: generate simulation data by specifying a certain type of missing mechanism.
#'
#' @docType package
#' @name missoNet-package
#' @useDynLib missoNet, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL
#> NULL