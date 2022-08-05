#' Fit missoNet models with a series of \eqn{\lambda} values
#'
#' @param X Numeric predictor matrix (\eqn{n\times p}{n x p}): columns correspond to predictor variables and rows correspond to samples. Missing values are not allowed. Do not include a column of ones.
#' @param Y Numeric response matrix (\eqn{n\times q}{n x q}): columns correspond to response variables and rows correspond to samples. Missing values should be coded as `\code{NA}`s or `\code{NaN}`s.
#' @param lambda.Beta A scalar or a numeric vector: an user-supplied sequence of non-negative \eqn{\lambda}{lambda} value(s) for regularizing the coefficient matrix \eqn{\mathbf{B}}{`\code{Beta}`}. Note that it must have a one-to-one correspondence with `\code{lambda.Theta}`.
#' @param lambda.Theta A scalar or a numeric vector: an user-supplied sequence of non-negative \eqn{\lambda}{lambda} value(s) for regularizing the precision matrix \eqn{\mathbf{\Theta}}{`\code{Theta}`}. Note that it must have a one-to-one correspondence with `\code{lambda.Beta}`.
#' @param rho (Optional) A scalar or a numeric vector of length \eqn{q}{q}: an user-supplied missing probability for response variables. Default is `\code{rho = NULL}` and the program will compute the empirical missing rates for columns of `\code{Y}` and use them as the working missing probability.
#' @param Beta.maxit The maximum number of iterations of the FISTA algorithm. Default is `\code{Beta.maxit = 10000}`.
#' @param Beta.thr The convergence threshold for updating \eqn{\mathbf{B}}{`\code{Beta}`}; default is `\code{Beta.thr = 1.0E-6}`. Iterations stop when absolute parameter change is less than `\code{Beta.thr * sum(abs(Beta))}`.
#' @param eta The backtracking line search shrinkage factor, default is `\code{eta = 0.8}`. Most users can use the default value, some experienced users may want to adjust `\code{eta}` according to the dataset's properties for a faster \eqn{\mathbf{B}}{`\code{Beta}`} convergence. Note that `\code{eta}` must be in (0, 1).
#' @param Theta.maxit The maximum number of iterations of the \code{glasso} algorithm. Default is `\code{Theta.maxit = 10000}`.
#' @param Theta.thr The convergence threshold for updating \eqn{\mathbf{\Theta}}{`\code{Theta}`}; default is `\code{Theta.thr = 1.0E-6}`. Iterations stop when average absolute parameter change is less than `\code{Theta.thr * ave(abs(offdiag(S)))}`, where `S` denotes the working empirical covariance matrix.
#' @param eps A numeric tolerance level for the L1 projection of the empirical covariance matrix; default is `\code{eps = 1.0E-8}`. The empirical covariance matrix will be projected onto a L1 ball to have `\code{min(eigen(S)$value) == eps}`, if any of the eigenvalues is less than the specified tolerance. Most users can use the default value.
#' @param penalize.diagonal Logical: should the diagonal of \eqn{\mathbf{\Theta}}{`\code{Theta}`} be penalized? The default depends on the sample size \eqn{n}{n} relative to the number of predictors and responses. If \eqn{n > \text{max}(p, q)}{n > max(p, q)}, the default is `TRUE`, otherwise it is set to `FALSE`.
#' @param diag.penalty.factor Numeric: a separate penalty factor for the diagonal entries of \eqn{\mathbf{\Theta}}{`\code{Theta}`} when `\code{penalize.diagonal = TRUE}`. \eqn{\lambda_\Theta}{`\code{lambda.Theta}`} is multiplied by this number to allow a differential shrinkage of the diagonal. The default is `\code{NULL}` and the program can compute it based on an initial estimate of \eqn{\mathbf{\Theta}}{`\code{Theta}`}. Could be `0` for no shrinkage (equivalent to `\code{penalize.diagonal = FALSE}`).
#' @param standardize Logical: should the columns of `\code{X}` be standardized so each has unit length? The default is `\code{TRUE}`. The estimated parameters will always be returned on the original scale. If `\code{X}` has been standardized prior to fitting the model, you might not wish to standardize.
#' @param standardize.response Logical: should the columns of `\code{Y}` be standardized so each has unit length? The default is `\code{TRUE}`. The estimated parameters will be returned on the original scale. If `\code{Y}` has been standardized prior to fitting the model, you might not wish to standardize.
#' @param fit.relax Logical: the default is `\code{FALSE}`. If `\code{TRUE}`, the program will re-estimate the edges (off-diagonal elements) in the active set of \eqn{\mathbf{\Theta}}{`\code{Theta}`} without penalization (\eqn{\lambda_\Theta=0}{`\code{lambda.Theta = 0}`}), which could be useful for further analyses of conditional inter-dependencies. WARNING: there may be convergence issues if the empirical covariance matrix is not of full rank (e.g., \eqn{n < q)}{n < q}).
#' @param parallel Logical: the default is `\code{FALSE}`. If `\code{TRUE}`, the program uses parallel clusters to fit models at each `\code{lambda}` value.
#' @param cpus Number of cores for parallelization. Only needed when `\code{parallel = TRUE}`.
#' @param verbose Value of `0`, `1` or `2`. `verbose = 0` -- silent; `verbose = 1` -- limited tracing; `verbose = 2` -- detailed tracing. Limited tracing if `\code{parallel = TRUE}`.
#'
#' @return This function returns a `\code{list}` of estimates for each \eqn{\lambda}{lambda} pair (\eqn{\lambda_B, \lambda_\Theta}{`\code{lambda.Beta}`, `\code{lambda.Theta}`}):
#' \item{\code{est.list}}{A `\code{list}` vector with the same length as \eqn{\lambda_B}{`\code{lambda.Beta}`}(\eqn{\lambda_\Theta}{`\code{lambda.Theta}`}). Each named `\code{list}` contains the following components:
#'   \itemize{
#'     \item \code{Beta}: the penalized estimate of the regression coefficient matrix (\eqn{p\times q}{p x q}).
#'     \item \code{Theta}: the penalized estimate of the precision matrix (\eqn{q\times q}{q x q}).
#'     \item \code{mu}: a vector of length \eqn{q}{q} storing the estimated intercept.
#'     \item \code{lambda.Beta}: the exact \eqn{\lambda_B}{`\code{lambda.Beta}`} value used to fit the model.
#'     \item \code{lambda.Theta}: the exact \eqn{\lambda_\Theta}{`\code{lambda.Theta}`} value used to fit the model.
#'     \item \code{relax.net}: a relaxed estimate of the conditional network structure (\eqn{q\times q}{q x q}) if `\code{fit.relax = TRUE}`.
#'   }
#' }
#' \item{\code{rho}}{A vector of length \eqn{q}{q} storing the working missing probability.}
#' \item{\code{penalize.diagonal}}{Logical: whether the diagonal of \eqn{\mathbf{\Theta}}{`\code{Theta}`} was penalized.}
#' \item{\code{diag.penalty.factor}}{The additional penalty factor for the diagonal of \eqn{\mathbf{\Theta}}{`\code{Theta}`} when `\code{penalize.diagonal = TRUE}`.}
#' @export
#'
#' @examples

missoNet <- function(X, Y, lambda.Beta, lambda.Theta, rho = NULL,
                     Beta.maxit = 1e4, Beta.thr = 1e-06, eta = 0.8,
                     Theta.maxit = 1e4, Theta.thr = 1e-06, eps = 1e-08,
                     penalize.diagonal = NULL, diag.penalty.factor = NULL,
                     standardize = TRUE, standardize.response = TRUE,
                     fit.relax = FALSE, parallel = FALSE, cpus = 2, verbose = 1) {
  if (length(lambda.Beta) != length(lambda.Theta)) {
    stop("`lambda.Beta` and `lambda.Theta` should be equal in length and have a one-to-one correspondence.")
  }
  
  if (verbose > 0) { cat("\n========================= missoNet ========================\n
- Parameter initialization...\n\n") }
  n <- nrow(X)
  kfold <- ifelse(n >= 100, 10, 5)
  set.seed(123)  # for reproducibility
  ind <- sample(n, replace = FALSE)
  foldid <- unlist(lapply(1:kfold, function(x) { rep(x, length((1 + floor((x - 1) * n/kfold)):floor(x * n/kfold))) }))
  init.obj <- InitParams(X = X[ind, ], Y = Y[ind, ], rho = rho, kfold = kfold, foldid = foldid,
                         Theta.maxit = Theta.maxit, Theta.thr = Theta.thr, eps = eps,
                         penalize.diagonal = penalize.diagonal, diag.pf = diag.penalty.factor,
                         standardize = standardize, standardize.response = standardize.response)
  
  if (verbose > 0) { cat("-----------------------------------------------------------\n
- Fittig with the `lambda` value(s)...\n\n") }
  if (!parallel) {
    fit_list <- vector("list", length(lambda.Theta))
    names(fit_list) <- paste0("pair", 1:length(lambda.Theta), ": lamB=", sprintf("%.2f", lambda.Beta), " lamTh=", sprintf("%.2f", lambda.Theta))
    if (verbose == 1) { pb <- txtProgressBar(min = 0, max = length(lambda.Theta), style = 3, width = 50, char = "=") }
    for (i in 1:length(lambda.Theta)) {
      fit_list[[i]] <- fitWrapper(X = X, Y = Y, lambda.Theta = lambda.Theta[i], lambda.Beta = lambda.Beta[i],
                                  Beta.maxit = Beta.maxit, Beta.thr = Beta.thr, eta = eta, Theta.maxit = Theta.maxit, Theta.thr = Theta.thr, eps = eps,
                                  verbose = verbose, fit.relax = fit.relax, init.obj = init.obj)
      if (verbose == 1) { setTxtProgressBar(pb, i) }
    }
    if (verbose == 1) { close(pb); cat("\n") }
  } else {
    if (verbose > 0) {
      snowfall::sfInit(parallel = TRUE, cpus = cpus)
    } else { suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = cpus)) }
    
    fit_list <- snowfall::sfClusterApplyLB(1:length(lambda.Theta), function(i) {
      fitWrapper(X = X, Y = Y, lambda.Theta = lambda.Theta[i], lambda.Beta = lambda.Beta[i],
                 Beta.maxit = Beta.maxit, Beta.thr = Beta.thr, eta = eta, Theta.maxit = Theta.maxit, Theta.thr = Theta.thr, eps = eps,
                 verbose = 0, fit.relax = fit.relax, init.obj = init.obj)
    })
    
    if (verbose > 0) {
      snowfall::sfStop()
    } else { suppressMessages(snowfall::sfStop()) }
    
    names(fit_list) <- paste0("pair", 1:length(lambda.Theta), ": lamB=", sprintf("%.2f", lambda.Beta), " lamTh=", sprintf("%.2f", lambda.Theta))
  }
  if (verbose > 0) { cat("========================= FINISHED ========================\n\n") }
  
  return(list(est.list = fit_list, rho = init.obj$rho.vec,
              penalize.diagonal = init.obj$penalize.diagonal, diag.penalty.factor = init.obj$diag.pf))
}


fitWrapper <- function(X, Y, lambda.Theta, lambda.Beta,
                       Beta.maxit, Beta.thr, eta, Theta.maxit, Theta.thr, eps,
                       verbose, fit.relax, init.obj) {
  fit <- update.missoNet(X = X, Y = Y, lamTh = lambda.Theta, lamB = lambda.Beta,
                         Beta.maxit = Beta.maxit, Beta.thr = Beta.thr,
                         Theta.maxit = Theta.maxit, Theta.thr = Theta.thr,
                         verbose = verbose, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                         info = NULL, info.update = NULL, init.obj = init.obj, under.cv = FALSE)
  fit$Beta <- sweep(fit$Beta/init.obj$sdx, 2, init.obj$sdy, `*`)    ## convert back to the original scale
  fit$mu <- as.numeric(init.obj$my - crossprod(fit$Beta, init.obj$mx))
  fit$lambda.Beta <- lambda.Beta
  fit$lambda.Theta <- lambda.Theta
  relax.net <- NULL
  if (fit.relax) {
    relax.net <- relax.glasso(X = X, Y = Y, init.obj = init.obj, est = fit, eps = eps,
                              Theta.thr = Theta.thr, Theta.maxit = Theta.maxit)
  }
  fit$relax.net <- relax.net
  
  return(fit)
}

