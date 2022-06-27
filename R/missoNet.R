#' Fit missoNet with a given pair of regularization parameters
#'
#' @param X Numeric predictor matrix (n by p): columns correspond to predictor variables and rows correspond to samples. Missing values are not allowed. Do not include a column of ones.
#' @param Y Numeric response matrix (n by q): columns correspond to response variables and rows correspond to samples. Missing values should be coded as \code{NA} or \code{NaN}.
#' @param lambda.Beta Numeric: a user supplied positive value for regularizing \code{Beta}.
#' @param lambda.Theta Numeric: a user supplied positive value for regularizing \code{Theta}.
#' @param rho (Optional) A scalar or a numeric vector of length q: a user supplied missing probability for response variables. Default is \code{rho = NULL} and the program will compute the empirical missing rates for columns of \code{Y} and use them as the missing probability.
#' @param Beta.maxit The maximum number of internal iterations allowed for updating \code{Beta}. Default is \code{Beta.maxit = 1e4}.
#' @param Beta.thr The convergence threshold for updating \code{Beta}; default is \code{Beta.thr = 1e-5}. Iterations stop when absolute parameter change is less than \code{Beta.thr * sum(abs(Beta))}.
#' @param eta Backtracking line search shrinkage factor, default is \code{eta = 0.8}. In some cases you may want to pick another \code{eta} for a faster \code{Beta} convergence based on your dataset. Note that \code{eta} must be positive and smaller than 1.
#' @param Theta.maxit The maximum number of internal iterations allowed for updating \code{Theta}. Default is \code{Theta.maxit = 1e4}.
#' @param Theta.thr The convergence threshold for updating \code{Theta}; default is \code{Theta.thr = 1e-5}. Iterations stop when average absolute parameter change is less than \code{Theta.thr * ave(abs(offdiag(Sigma)))}.
#' @param eps A numeric tolerance level for L1 projection; default is \code{eps = 1e-8}. If any of the eigenvalues is less than the given tolerance, the unbiased estimate of covariance is projected onto L1 ball to have \code{min(eigen(Sigma)$value) == eps}.
#' @param diag.penalty.factor Numeric: a separate penalty factor for the diagonal entries of \code{Theta}. This is a number that multiplies \code{lambda.Theta} to allow differential shrinkage. Default is \code{NULL} and the program computes it based on an initial estimate of \code{Theta}. Can be 0 for no shrinkage. Only needed when n <= p.
#' @param standardize Logical: should the columns of \code{X} be standardized so each has unit length and zero average; default is \code{TRUE}. The parameter estimates will be returned on the original scale. If \code{X} has been standardized prior to fitting the model, you might not wish to standardize.
#' @param standardize.response Logical: should the columns of \code{Y} be standardized so each has unit length and zero average; default is \code{TRUE}. The parameter estimates will be returned on the original scale. If \code{Y} has been standardized prior to fitting the model, you might not wish to standardize.
#' @param fit.relax Logical: default is \code{FALSE}. If \code{TRUE}, the program will re-estimate a relaxed graph (\code{Theta}) without penalization, which could be useful for network analysis. WARNING: there may be convergence issues if the residual covariance matrix is not of full rank.
#' @param verbose Value of 0, 1 or 2. 0 -- silent; 1 -- limited tracing; 2 -- detailed tracing.
#'
#' @return This function returns a \code{list} of estimates and working parameters:
#' \itemize{
#' \item \code{est}: estimates using the designated regularization parameters.
#' \item \code{relax.graph}: an estimate of network without penalization if \code{fit.relax = TRUE}.
#' \item \code{rho}: a vector of length q: working missing probability.
#' \item \code{penalize.diagonal}: logical: whether the diagonal of \code{Theta} was penalized.
#' \item \code{diag.penalty.factor}: penalty factor for the diagonal when \code{penalize.diagonal} returned as \code{TRUE}.
#' }
#' @export
#'
#' @examples

missoNet <- function(X, Y, lambda.Beta, lambda.Theta, rho = NULL,
                     Beta.maxit = 1e4, Beta.thr = 1e-05, eta = 0.8, Theta.maxit = 1e4, Theta.thr = 1e-05, eps = 1e-08,
                     diag.penalty.factor = NULL, standardize = TRUE, standardize.response = TRUE, fit.relax = FALSE, verbose = 1) {
  if (verbose > 0) {
    cat("\nInitializing necessary parameters...\n\n")
  }
  
  n <- nrow(X)
  kfold <- ifelse(n >= 100, 10, 5)
  set.seed(NULL)
  ind <- sample(n, replace = FALSE)
  foldid <- unlist(lapply(1:kfold, function(x) {
    rep(x, length((1 + floor((x - 1) * n/kfold)):floor(x * n/kfold)))
  }))
  
  init.obj <- InitParams(X = X[ind, ], Y = Y[ind, ], rho = rho, kfold = kfold, foldid = foldid,
                         Theta.maxit = Theta.maxit, Theta.thr = Theta.thr, eps = eps, diag.pf = diag.penalty.factor,
                         standardize = standardize, standardize.response = standardize.response)
  
  if (verbose > 0) {
    cat("Fittig with the given lambda pair...\n\n")
  }
  fit <- update.missoNet(X = X, Y = Y, lamTh = lambda.Theta, lamB = lambda.Beta,
                         Beta.maxit = Beta.maxit, Beta.thr = Beta.thr,
                         Theta.maxit = Theta.maxit, Theta.thr = Theta.thr,
                         verbose = verbose, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                         info = NULL, info.update = NULL, init.obj = init.obj, under.cv = FALSE)
  fit$lambda.Beta <- lambda.Beta
  fit$lambda.Theta <- lambda.Theta
  fit$Beta <- sweep(fit$Beta/init.obj$sdx, 2, init.obj$sdy, `*`)    ## convert back to the original scale
  fit$mu <- as.numeric(init.obj$my - crossprod(fit$Beta, init.obj$mx))
  if (verbose > 0) {
    cat("Done.\n\n")
  }
  
  relax.graph <- NULL
  if (fit.relax) {
    relax.graph <- relax.glasso(X = X, Y = Y, init.obj = init.obj, est = fit, eps = eps,
                                Theta.thr = Theta.thr, Theta.maxit = Theta.maxit)
  }
  
  return(list(est = fit, relax.graph = relax.graph, rho = init.obj$rho.vec, 
              penalize.diagonal = init.obj$penalize.diagonal, diag.penalty.factor = init.obj$diag.pf))
}


