#' Make predictions from a \code{cv.missoNet} object
#'
#' @param object Fitted \code{cv.missoNet} object.
#' @param newx Predictor matrix of new values at which predictions are to be made. \code{newx} should be in the same scale as the training data. Missing values are not allowed. Do not include a column of ones.
#' @param s The penalty parameter lambda at which the regression coefficients for predictions are extracted. Default is \code{s = "lambda.min"} stored in the \code{cv.missoNet} object. Alternatively \code{s = "lambda.1se"} can be used. 
#' @param ... Not used. Other arguments for predicting.
#' 
#' @return Predicted response matrix \code{newy}.
#' 
#' @method predict cv.missoNet
#' @export
#'
#' @examples

predict.cv.missoNet <- function(object, newx = NULL, s = "lambda.min", ...) {
  if (is.null(newx)) {
    stop("please supply a predictor matrix (n by p) in the same scale as the training data.")
  }
  n <- nrow(newx)
  if (s == "lambda.min") {
    return(rep(1, n) %*% t(object$est.min$mu) + newx %*% object$est.min$Beta)
  } else if (s == "lambda.1se") {
    if (is.null(object$est.1se)) {
      stop("please train the model by setting 'fit.1se = TRUE'.")
    } else {
      return(rep(1, n) %*% t(object$est.1se$mu) + newx %*% object$est.1se$Beta)
    }
  } else {
    stop("invalid form for s, s should be either 'lambda.min' or 'lambda.1se'.")
  }
}