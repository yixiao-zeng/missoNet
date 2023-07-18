#' Make predictions from a cv.missoNet object
#'
#' S3 method for making predictions of response values from a fitted \code{'cv.missoNet'} object.
#'
#' @param object A fitted \code{'cv.missoNet'} object.
#' @param newx A predictor matrix of new values at which predictions are to be made. The columns of \code{'newx'} should have the same standardization flags as the original input for training the model. Missing values are not allowed. \code{'newx'} should not include a column of ones for an intercept.
#' @param s Character string, the regularization parameter pair \eqn{\lambda} = (\eqn{\lambda_B}, \eqn{\lambda_\Theta}) at which the coefficients are extracted for making predictions. It supports three special strings, named "\code{lambda.min}" (default), "\code{lambda.1se.Beta}" and "\code{lambda.1se.Theta}". 
#' @param ... Not used. Other arguments for predicting.
#' 
#' @return The matrix of predicted values: \code{'newy = mu_hat + newx \%*\% Beta_hat'}.
#' @method predict cv.missoNet
#' @export
#' 
#' @author Yixiao Zeng \email{yixiao.zeng@@mail.mcgill.ca}, Celia M.T. Greenwood and Archer Yi Yang.
#' 
#' @examples
#' ## Simulate a dataset.
#' set.seed(123)  # reproducibility
#' sim.dat <- generateData(n = 300, p = 10, q = 10, rho = 0.1, missing.type = "MCAR")
#' tr <- 1:240  # training set indices
#' tst <- 241:300  # test set indices
#' 
#' \donttest{
#' ## Perform a five-fold cross-validation on the training set.
#' cvfit <- cv.missoNet(X = sim.dat$X[tr, ], Y = sim.dat$Z[tr, ], kfold = 5,
#'                      fit.1se = TRUE, permute = TRUE, with.seed = 486)
#' 
#' 
#' ## Make predictions of response values on the test set.
#' newy1 <- predict(cvfit, newx = sim.dat$X[tst, ], s = "lambda.min")
#' newy2 <- predict(cvfit, newx = sim.dat$X[tst, ], s = "lambda.1se.Beta")  # 'fit.1se' = TRUE
#' newy3 <- predict(cvfit, newx = sim.dat$X[tst, ], s = "lambda.1se.Theta")  # 'fit.1se' = TRUE
#' }

predict.cv.missoNet <- function(object, newx = NULL, s = "lambda.min", ...) {
  if (is.null(newx)) {
    stop("\nPlease supply a predictor matrix (n* x p) for `newx`.\n")
  }
  n <- nrow(newx)
  if (s == "lambda.min") {
    return(rep(1, n) %*% t(object$est.min$mu) + newx %*% object$est.min$Beta)
  } else if (s == "lambda.1se.Beta") {
    if (is.null(object$est.1se.B)) {
      stop("\n`lambda.1se.Beta` not found. Please make sure that `fit.1se = TRUE` and `est.1se.B` is not `NULL`.\n")
    } else {
      return(rep(1, n) %*% t(object$est.1se.B$mu) + newx %*% object$est.1se.B$Beta)
    }
  } else if (s == "lambda.1se.Theta") {
    if (is.null(object$est.1se.Tht)) {
      stop("\n`lambda.1se.Theta` not found. Please make sure that `fit.1se = TRUE` and `est.1se.Tht` is not `NULL`.\n")
    } else {
      return(rep(1, n) %*% t(object$est.1se.Tht$mu) + newx %*% object$est.1se.Tht$Beta)
    }
  } else {
    stop('\nInvalid input for `s`, `s` should be one of `"lambda.min"`, `"lambda.1se.Beta"` or `"lambda.1se.Theta"`.\n')
  }
}

