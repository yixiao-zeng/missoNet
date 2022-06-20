#' Make predictions from a cv.missoNet object.
#'
#' @param object Fitted cv.missoNet object.
#' @param newx Matrix of new values for X at which predictions are to be made.
#' @param s The penalty parameter lambda at which predictions are required. Default is \code{s = "lambda.min"} stored in the CV object. Alternatively \code{s = "lambda.1se"} can be used. 
#' @param ... Not used. Other arguments for predicting.
#' 
#' @return Predicted response matrix code{newY}.
#' 
#' @method predict cv.missoNet
#' @export
#'
#' @examples

predict.cv.missoNet <- function(object, newx = NULL, s = "lambda.min", ...) {
  if (is.null(newx)) {
    stop("please supply a predictor matrix (n by p) having the same scale as the training data.")
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