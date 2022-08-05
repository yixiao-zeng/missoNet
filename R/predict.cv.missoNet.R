#' Make predictions from a \code{cv.missoNet} object
#'
#' S3 method for making predictions from a fitted `\code{cv.missoNet}` object.
#'
#' @param object Fitted `\code{cv.missoNet}` object.
#' @param newx A predictor matrix of new values at which predictions are to be made. `\code{newx}` should have the same scale as the training data. Missing values are not allowed and do not include a column of ones.
#' @param s Character, the penalty parameter `lambda` at which the regression coefficients for predictions are extracted. It supports three special values, named `\code{"lambda.min"}` (default), `\code{"lambda.Beta.1se"}` and `\code{"lambda.Theta.1se"}`. 
#' @param ... Not used. Other arguments for predicting.
#' 
#' @return Predicted response matrix `\code{newy}`.
#' 
#' @method predict cv.missoNet
#' @export

predict.cv.missoNet <- function(object, newx = NULL, s = "lambda.min", ...) {
  if (is.null(newx)) {
    stop("\nPlease supply a predictor matrix (n by p) having the same scale as the training data.\n")
  }
  n <- nrow(newx)
  if (s == "lambda.min") {
    return(rep(1, n) %*% t(object$est.min$mu) + newx %*% object$est.min$Beta)
  } else if (s == "lambda.Beta.1se") {
    if (is.null(object$estB.1se)) {
      stop("\n`lambda.Beta.1se` not found. Please make sure `fit.1se = TRUE` and `estB.1se` is not `NULL`.\n")
    } else {
      return(rep(1, n) %*% t(object$estB.1se$mu) + newx %*% object$estB.1se$Beta)
    }
  } else if (s == "lambda.Theta.1se") {
    if (is.null(object$estTht.1se)) {
      stop("\n`lambda.Theta.1se` not found. Please make sure `fit.1se = TRUE` and `estTht.1se` is not `NULL`.\n")
    } else {
      return(rep(1, n) %*% t(object$estTht.1se$mu) + newx %*% object$estTht.1se$Beta)
    }
  } else {
    stop('\nInvalid input for `s`, `s` should be one of `"lambda.min"`, `"lambda.Beta.1se"` or `"lambda.Theta.1se"`.\n')
  }
}

