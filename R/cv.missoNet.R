#' Cross-validation for missoNet
#'
#' @param X Numeric predictor matrix (\eqn{n\times p}{n x p}): columns correspond to predictor variables and rows correspond to samples. Missing values are not allowed. Do not include a column of ones.
#' @param Y Numeric response matrix (\eqn{n\times q}{n x q}): columns correspond to response variables and rows correspond to samples. Missing values should be coded as `\code{NA}`s or `\code{NaN}`s.
#' @param kfold Number of folds -- the default is `\code{5}`.
#' @param rho (Optional) A scalar or a numeric vector of length \eqn{q}{q}: an user-supplied missing probability for response variables. Default is `\code{rho = NULL}` and the program will compute the empirical missing rates for columns of `\code{Y}` and use them as the working missing probability.
#' @param lambda.Beta (Optional) Numeric vector: an user-supplied sequence of non-negative \eqn{\lambda}{`\code{lambda}`} values for regularizing \eqn{\mathbf{B}}{`\code{Beta}`} from which the CV procedure searches. Default is `\code{lambda.Beta = NULL}` and the program computes an appropriate \eqn{\lambda_B}{`\code{lambda.Beta}`} sequence based on `\code{n.lamBeta}` and `\code{lamBeta.min.ratio}`, supplying a vector overrides this. The sequence supplied will be automatically arranged in a descending order internally. 
#' @param lambda.Theta (Optional) Numeric vector: an user-supplied sequence of non-negative \eqn{\lambda}{`\code{lambda}`} values for regularizing \eqn{\mathbf{\Theta}}{`\code{Theta}`} from which the CV procedure searches. Default is `\code{lambda.Theta = NULL}` and the program computes an appropriate \eqn{\lambda_\Theta}{`\code{lambda.Theta}`} sequence based on `\code{n.lamTheta}` and `\code{lamTheta.min.ratio}`, supplying a vector overrides this. The sequence supplied will be automatically arranged in a descending order internally. 
#' @param lamBeta.min.ratio The smallest value of \eqn{\lambda_B}{`\code{lambda.Beta}`} will be the data-derived `\code{lambda.Beta.max}` multiplied by `\code{lamBeta.min.ratio}`. The default depends on the sample size \eqn{n}{n} relative to the number of predictors \eqn{p}{p}. If \eqn{n > p}{n > p}, the default is `1.0E-4`, otherwise it is `1.0E-2`. A very small value of `\code{lamBeta.min.ratio}` may significantly increase the runtime and lead to a saturated fit in the \eqn{n \leq p}{n <= p} case. Only needed when `\code{lambda.Beta = NULL}`.
#' @param lamTheta.min.ratio The smallest value of \eqn{\lambda_\Theta}{`\code{lambda.Theta}`} will be the data-derived `\code{lambda.Theta.max}` multiplied by `\code{lamTheta.min.ratio}`. The default depends on the sample size \eqn{n}{n} relative to the number of responses \eqn{q}{q}. If \eqn{n > q}{n > q}, the default is `1.0E-4`, otherwise it is `1.0E-2`. A very small value of `\code{lamTheta.min.ratio}` may significantly increase the runtime and lead to a saturated fit in the \eqn{n \leq q}{n <= q} case. Only needed when `\code{lambda.Theta = NULL}`.
#' @param n.lamBeta The number of \eqn{\lambda_B}{`\code{lambda.Beta}`} values. If \eqn{n > p}{n > p}, the default is `40`, otherwise it is `20`. Avoid supplying a huge number since the program will fit `n.lamBeta * n.lamTheta` models in total for each fold of CV, typically we suggest `n.lamBeta = -log10(lamBeta.min.ratio) * c`, where `c` \eqn{\in} [10, 20]. Only needed when `\code{lambda.Beta = NULL}`.
#' @param n.lamTheta The number of \eqn{\lambda_\Theta}{`\code{lambda.Theta}`} values. If \eqn{n > q}{n > q}, the default is `40`, otherwise it is `20`. Avoid supplying a huge number since the program will fit `n.lamBeta * n.lamTheta` models in total for each fold of CV, typically we suggest `n.lamTheta = -log10(lamTheta.min.ratio) * c`, where `c` \eqn{\in} [10, 20]. Only needed when `\code{lambda.Theta = NULL}`.
#' @param lamBeta.scale.factor A positive multiplication factor for scaling the entire \eqn{\lambda_B}{`\code{lambda.Beta}`} sequence; the default is `1`. A typical usage scenario is when the optimal selection of \eqn{\lambda_B}{`\code{lambda.Beta}`} approaches the boundaries of the search range. Only needed when `\code{lambda.Beta = NULL}`.
#' @param lamTheta.scale.factor A positive multiplication factor for scaling the entire \eqn{\lambda_\Theta}{`\code{lambda.Theta}`} sequence; the default is `1`. A typical usage scenario is when the optimal selection of \eqn{\lambda_\Theta}{`\code{lambda.Theta}`} approaches the boundaries of the search range. Only needed when `\code{lambda.Theta = NULL}`.
#' @param Beta.maxit The maximum number of iterations of the FISTA algorithm. Default is `\code{Beta.maxit = 1000}`.
#' @param Beta.thr The convergence threshold for updating \eqn{\mathbf{B}}{`\code{Beta}`}; default is `\code{Beta.thr = 1.0E-4}`. Iterations stop when absolute parameter change is less than `\code{Beta.thr * sum(abs(Beta))}`.
#' @param eta The backtracking line search shrinkage factor, default is `\code{eta = 0.8}`. Most users can use the default value, some experienced users may want to adjust `\code{eta}` according to the dataset's properties for a faster \eqn{\mathbf{B}}{`\code{Beta}`} convergence. Note that `\code{eta}` must be in (0, 1).
#' @param Theta.maxit The maximum number of iterations of the \code{glasso} algorithm. Default is `\code{Theta.maxit = 1000}`.
#' @param Theta.thr The convergence threshold for updating \eqn{\mathbf{\Theta}}{`\code{Theta}`}; default is `\code{Theta.thr = 1.0E-4}`. Iterations stop when average absolute parameter change is less than `\code{Theta.thr * ave(abs(offdiag(S)))}`, where `S` denotes the working empirical covariance matrix.
#' @param eps A numeric tolerance level for the L1 projection of the empirical covariance matrix; default is `\code{eps = 1.0E-8}`. The empirical covariance matrix will be projected onto a L1 ball to have `\code{min(eigen(S)$value) == eps}`, if any of the eigenvalues is less than the specified tolerance. Most users can use the default value.
#' @param penalize.diagonal Logical: should the diagonal of \eqn{\mathbf{\Theta}}{`\code{Theta}`} be penalized? The default depends on the sample size \eqn{n}{n} relative to the number of predictors and responses. If \eqn{n > \text{max}(p, q)}{n > max(p, q)}, the default is `TRUE`, otherwise it is set to `FALSE`.
#' @param diag.penalty.factor Numeric: a separate penalty factor for the diagonal entries of \eqn{\mathbf{\Theta}}{`\code{Theta}`} when `\code{penalize.diagonal = TRUE}`. \eqn{\lambda_\Theta}{`\code{lambda.Theta}`} is multiplied by this number to allow a differential shrinkage of the diagonal. The default is `\code{NULL}` and the program can compute it based on an initial estimate of \eqn{\mathbf{\Theta}}{`\code{Theta}`}. Could be `0` for no shrinkage (equivalent to `\code{penalize.diagonal = FALSE}`).
#' @param standardize Logical: should the columns of `\code{X}` be standardized so each has unit length? The default is `\code{TRUE}`. The estimated parameters will always be returned on the original scale. If `\code{X}` has been standardized prior to fitting the model, you might not wish to standardize.
#' @param standardize.response Logical: should the columns of `\code{Y}` be standardized so each has unit length? The default is `\code{TRUE}`. The estimated parameters will be returned on the original scale. If `\code{Y}` has been standardized prior to fitting the model, you might not wish to standardize.
#' @param fit.1se Logical: should the model be refitted at the largest \eqn{\lambda_B}{`\code{lambda.Beta}`} and \eqn{\lambda_\Theta}{`\code{lambda.Theta}`} according to the one-standard-error rule? The default is `\code{FALSE}`.
#' @param fit.relax Logical: the default is `\code{FALSE}`. If `\code{TRUE}`, the program will re-estimate the edges (off-diagonal elements) in the active set of \eqn{\mathbf{\Theta}}{`\code{Theta}`} without penalization (\eqn{\lambda_\Theta=0}{`\code{lambda.Theta = 0}`}), which could be useful for further analyses of conditional inter-dependencies. WARNING: there may be convergence issues if the empirical covariance matrix is not of full rank (e.g., \eqn{n < q)}{n < q}).
#' @param permute Logical: should the subject indices for CV be permuted? The default is `\code{TRUE}`.
#' @param with.seed A random seed for permutation.
#' @param parallel Logical: the default is `\code{FALSE}`. If `\code{TRUE}`, the program uses parallel clusters to fit each fold of CV.
#' @param cpus Number of cores for parallelization. Only needed when `\code{parallel = TRUE}`.
#' @param verbose Value of `0`, `1` or `2`. `verbose = 0` -- silent; `verbose = 1` -- limited tracing; `verbose = 2` -- detailed tracing. Limited tracing if `\code{parallel = TRUE}`.
#'
#' @return This function returns a `\code{cv.missoNet}` object containing a named `\code{list}` with all the ingredients of the cross-validated fit:
#' \item{\code{est.min}}{A list of parameters estimated at `\code{lambda.min}` that gives the smallest mean CV error. It consists of the following components:
#'   \itemize{
#'       \item \code{Beta}: the penalized estimate of the regression coefficient matrix (\eqn{p\times q}{p x q}).
#'       \item \code{Theta}: the penalized estimate of the precision matrix (\eqn{q\times q}{q x q}).
#'       \item \code{mu}: a vector of length \eqn{q}{q} storing the estimated intercept.
#'       \item \code{lambda.Beta}: the exact \eqn{\lambda_B}{`\code{lambda.Beta}`} value used to fit the model.
#'       \item \code{lambda.Theta}: the exact \eqn{\lambda_\Theta}{`\code{lambda.Theta}`} value used to fit the model.
#'       \item \code{relax.net}: a relaxed estimate of the conditional network structure (\eqn{q\times q}{q x q}) if `\code{fit.relax = TRUE}`.
#'   }
#' }
#' \item{\code{estB.1se}}{A list of parameters estimated at `\code{lambda.Beta.1se}` if `\code{fit.1se = TRUE}`. `\code{lambda.Beta.1se}` is the largest \eqn{\lambda_B}{`\code{lambda.Beta}`} at which the mean CV error is within one standard error of the minimum, by fixing \eqn{\lambda_\Theta}{`\code{lambda.Theta}`} at `\code{lambda.min}`. It consists of the same components as `\code{est.min}`.}
#' \item{\code{estTht.1se}}{A list of parameters estimated at `\code{lambda.Theta.1se}` if `\code{fit.1se = TRUE}`. `\code{lambda.Theta.1se}` is the largest \eqn{\lambda_\Theta}{`\code{lambda.Theta}`} at which the mean CV error is within one standard error of the minimum, by fixing \eqn{\lambda_B}{`\code{lambda.Beta}`} at `\code{lambda.min}`. It consists of the same components as `\code{est.min}`.}
#' \item{\code{rho}}{A vector of length \eqn{q}{q} storing the working missing probability.}
#' \item{\code{fold.index}}{The subject indices identifying which fold each observation is in.}
#' \item{\code{lambda.Beta.vec}}{A long vector of length `\code{n.lamBeta * n.lamTheta}` from which the CV procedure searches.}
#' \item{\code{lambda.Theta.vec}}{A long vector of length `\code{n.lamBeta * n.lamTheta}` from which the CV procedure searches.}
#' \item{\code{cvm}}{A long vector of standardized mean CV error that has a one to one correspondence with `\code{lambda.Beta.vec}` and `\code{lambda.Theta.vec}`.}
#' \item{\code{cvup}}{Upper cross-validated error.}
#' \item{\code{cvlo}}{Lower cross-validated error.}
#' \item{\code{penalize.diagonal}}{Logical: whether the diagonal of \eqn{\mathbf{\Theta}}{`\code{Theta}`} was penalized.}
#' \item{\code{diag.penalty.factor}}{The additional penalty factor for the diagonal of \eqn{\mathbf{\Theta}}{`\code{Theta}`} when `\code{penalize.diagonal = TRUE}`.}
#' @export
#' 
#' @examples 

cv.missoNet <- function(X, Y, kfold = 5, rho = NULL,
                        lambda.Beta = NULL, lambda.Theta = NULL,
                        lamBeta.min.ratio = NULL, lamTheta.min.ratio = NULL,
                        n.lamBeta = NULL, n.lamTheta = NULL,
                        lamBeta.scale.factor = 1, lamTheta.scale.factor = 1,
                        Beta.maxit = 1000, Beta.thr = 1e-04, eta = 0.8,
                        Theta.maxit = 1000, Theta.thr = 1e-04, eps = 1e-08,
                        penalize.diagonal = NULL, diag.penalty.factor = NULL,
                        standardize = TRUE, standardize.response = TRUE,
                        fit.1se = FALSE, fit.relax = FALSE,
                        permute = TRUE, with.seed = NULL,
                        parallel = FALSE, cpus = 2, verbose = 1) {
  if (verbose > 0) { cat("\n======================= cv.missoNet =======================\n
- Parameter initialization...\n\n") }
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  if (permute) {
    set.seed(with.seed)
    ind <- sample(n, replace = FALSE)
  } else { ind <- 1:n }
  foldid <- unlist(lapply(1:kfold, function(x) { rep(x, length((1 + floor((x - 1) * n/kfold)):floor(x * n/kfold))) }))
  names(ind) <- paste0("fold.", foldid)
  
  init.obj <- InitParams(X = X[ind, ], Y = Y[ind, ], rho = rho, kfold = kfold, foldid = foldid,
                         Theta.maxit = Theta.maxit, Theta.thr = Theta.thr, eps = eps, 
                         penalize.diagonal = penalize.diagonal, diag.pf = diag.penalty.factor,
                         standardize = standardize, standardize.response = standardize.response)
  
  lambda.obj <- InitLambda(lamB = lambda.Beta, lamTh = lambda.Theta, n.tr = floor(n * (kfold - 1)/kfold),
                           init.obj = init.obj, n.lamB = n.lamBeta, n.lamTh = n.lamTheta,
                           lamB.min.ratio = lamBeta.min.ratio, lamTh.min.ratio = lamTheta.min.ratio,
                           lamB.scale.factor = lamBeta.scale.factor, lamTh.scale.factor = lamTheta.scale.factor)
  lamTh.vec <- lambda.obj$lamTh.vec
  lamB.vec <- lambda.obj$lamB.vec
  
  ################################################################################
  # Cross-validation
  ################################################################################
  if (verbose > 0) { cat("--------------------- Cross-validation --------------------\n\n") }
  if (!parallel) {
    err <- matrix(0, kfold, length(lamTh.vec))
    for (k in 1:kfold) {
      if (verbose == 1) { cat(sprintf("Fold: %d/%d\n", k, kfold)) }
      
      foldind <- ind[(1 + floor((k - 1) * n/kfold)):floor(k * n/kfold)]
      X.tr <- X[-foldind, ]
      Y.tr <- Y[-foldind, ]
      X.va <- X[foldind, , drop = FALSE]
      Y.va <- Y[foldind, , drop = FALSE]
      n.tr <- nrow(Y.tr)
      
      if (is.null(rho)) {
        rho.vec <- apply(Y.tr, 2, function(x) {
          sum(is.na(x))/n.tr
        })
      } else { rho.vec <- init.obj$rho.vec }
      rho.mat.1 <- t(matrix(rep(1 - rho.vec, p), q, p))  ## pxq
      rho.mat.2 <- matrix(1 - rho.vec, q, 1) %*% matrix(1 - rho.vec, 1, q)
      diag(rho.mat.2) <- 1 - rho.vec  ## qxq
      
      mx.tr <- apply(X.tr, 2, mean)
      X.tr <- scale(X.tr, center = mx.tr, scale = init.obj$sdx)
      X.va <- scale(X.va, center = mx.tr, scale = init.obj$sdx)
      
      my.tr <- apply(Y.tr, 2, mean, na.rm = TRUE)
      Y.tr <- scale(Y.tr, center = my.tr, scale = init.obj$sdy)
      Y.va <- scale(Y.va, center = my.tr, scale = init.obj$sdy)
      Z.tr <- Y.tr
      Z.tr[is.na(Z.tr)] <- 0
      
      info <- NULL
      info$n <- n.tr
      info$q <- q
      info$penalize.diagonal <- init.obj$penalize.diagonal
      info$xtx <- crossprod(X.tr)
      info$til.xty <- crossprod(X.tr, Z.tr)/rho.mat.1
      
      info.update <- NULL
      info.update$B.init <- init.obj$B.init * init.obj$sdx   ## initialize B on the standardized scale
      Beta.thr.rescale <- Beta.thr * sum(abs(info.update$B.init))
      E.tr <- Y.tr - X.tr %*% info.update$B.init
      info.update$residual.cov <- getResidual(E = E.tr, n = n.tr, rho.mat = rho.mat.2, eps = eps)
      
      if (verbose == 1) { pb <- txtProgressBar(min = 0, max = length(lamTh.vec), style = 3, width = 50, char = "=") }
      for (i in 1:length(lamTh.vec)) {
        info.update$B.init <- update.missoNet(lamTh = lamTh.vec[i], lamB = lamB.vec[i],
                                              Beta.maxit = Beta.maxit, Beta.thr = Beta.thr.rescale,
                                              Theta.maxit = Theta.maxit, Theta.thr = Theta.thr,
                                              verbose = verbose, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                                              info = info, info.update = info.update, under.cv = TRUE)
        Beta.thr.rescale <- Beta.thr * sum(abs(info.update$B.init))
        E.tr <- Y.tr - X.tr %*% info.update$B.init
        info.update$residual.cov <- getResidual(E = E.tr, n = n.tr, rho.mat = rho.mat.2, eps = eps)
        
        E.va.sq <- (Y.va - X.va %*% info.update$B.init)^2
        err[k, i] <- mean(E.va.sq, na.rm = TRUE)
        if (verbose == 1) { setTxtProgressBar(pb, i) }
      }
      if (verbose == 1) { close(pb) }
    }
    if (verbose == 1) { cat("\n") }

  } else {
    if (verbose > 0) {
      snowfall::sfInit(parallel = TRUE, cpus = cpus)
    } else { suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = cpus)) }
    
    par.out <- snowfall::sfClusterApplyLB(seq(1, kfold, by = 1), function(k) {
      parWrapper(k = k, X = X, Y = Y, init.obj = init.obj, rho = rho, ind = ind, kfold = kfold, lamTh.vec = lamTh.vec, lamB.vec = lamB.vec,
                 Beta.maxit = Beta.maxit, Beta.thr = Beta.thr, Theta.maxit = Theta.maxit, Theta.thr = Theta.thr, eps = eps, eta = eta)
    })
    
    if (verbose > 0) {
      snowfall::sfStop()
    } else { suppressMessages(snowfall::sfStop()) }
    
    err <- do.call("rbind", par.out)
  }
  if (verbose > 0) { cat("-----------------------------------------------------------\n\n") }
  
  err.cv <- colSums(err)/kfold
  err.sd <- apply(err, 2, sd)/sqrt(kfold)
  err.up <- err.cv + err.sd
  err.low <- err.cv - err.sd
  
  cv.min <- which.min(err.cv)
  lamTh.min <- lamTh.vec[cv.min]
  lamB.min <- lamB.vec[cv.min]
  boundaryCheck(lambda.Theta = lambda.Theta, lambda.Beta = lambda.Beta, 
                lamTh.vec = lamTh.vec, lamB.vec = lamB.vec, 
                lamTh.min = lamTh.min, lamB.min = lamB.min, margin = 0.1)
  
  if (verbose > 0) { cat("- Fittig with `lambda.min`...\n\n") }
  out.min <- update.missoNet(X = X, Y = Y, lamTh = lamTh.min, lamB = lamB.min,
                             Beta.maxit = Beta.maxit * 10, Beta.thr = Beta.thr * 0.01,
                             Theta.maxit = Theta.maxit * 10, Theta.thr = Theta.thr * 0.01,
                             verbose = verbose, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                             info = NULL, info.update = NULL, init.obj = init.obj, under.cv = FALSE)
  out.min$Beta <- sweep(out.min$Beta/init.obj$sdx, 2, init.obj$sdy, `*`)    ## convert back to the original scale
  out.min$mu <- as.numeric(init.obj$my - crossprod(out.min$Beta, init.obj$mx))
  out.min$lambda.Beta <- lamB.min
  out.min$lambda.Theta <- lamTh.min
  relax.net <- NULL
  if (fit.relax) {
    if (verbose > 0) { cat("  -- Relaxed network inference...\n\n") }
    relax.net <- relax.glasso(X = X, Y = Y, init.obj = init.obj, est = out.min, eps = eps,
                              Theta.thr = Theta.thr * 0.01, Theta.maxit = Theta.maxit * 10)
  }
  out.min$relax.net <- relax.net
  
  outB.1se <- NULL
  outTh.1se <- NULL
  if (fit.1se) {
    if (verbose > 0) { 
      cat("------------------------------------\n
- Fiting with `lambda.1se`...\n
  -- `lambda.Beta.1se`\n\n") }
    new.lamB.vec <- lamB.vec[lamTh.vec == lamTh.min]
    new.err.cv <- err.cv[lamTh.vec == lamTh.min]
    lamB.1se <- max(new.lamB.vec[new.err.cv <= err.up[cv.min]])
    
    if (lamB.1se != lamB.min) {
      outB.1se <- update.missoNet(X = X, Y = Y, lamTh = lamTh.min, lamB = lamB.1se,
                                  Beta.maxit = Beta.maxit * 10, Beta.thr = Beta.thr * 0.01,
                                  Theta.maxit = Theta.maxit * 10, Theta.thr = Theta.thr * 0.01,
                                  verbose = verbose, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                                  info = NULL, info.update = NULL, init.obj = init.obj, under.cv = FALSE)
      outB.1se$Beta <- sweep(outB.1se$Beta/init.obj$sdx, 2, init.obj$sdy, `*`)
      outB.1se$mu <- as.numeric(init.obj$my - crossprod(outB.1se$Beta, init.obj$mx))
      outB.1se$lambda.Beta <- lamB.1se
      outB.1se$lambda.Theta <- lamTh.min
      relax.net <- NULL
      if (fit.relax) {
        if (verbose > 0) { cat("    --- Relaxed network inference...\n\n") }
        relax.net <- relax.glasso(X = X, Y = Y, init.obj = init.obj, est = outB.1se, eps = eps,
                                  Theta.thr = Theta.thr * 0.01, Theta.maxit = Theta.maxit * 10)
      }
      outB.1se$relax.net <- relax.net
    } else {
      warning("\n`lambda.Beta.1se` = `lambda.Beta.min`, please provide a finer grid for `lambda.Beta` by increasing the number of values.\n")
    }
    
    if (verbose > 0) { cat("  -- `lambda.Theta.1se`\n\n") }
    new.lamTh.vec <- lamTh.vec[lamB.vec == lamB.min]
    new.err.cv <- err.cv[lamB.vec == lamB.min]
    lamTh.1se <- max(new.lamTh.vec[new.err.cv <= err.up[cv.min]])
    
    if (lamTh.1se != lamTh.min) {
      outTh.1se <- update.missoNet(X = X, Y = Y, lamTh = lamTh.1se, lamB = lamB.min,
                                  Beta.maxit = Beta.maxit * 10, Beta.thr = Beta.thr * 0.01,
                                  Theta.maxit = Theta.maxit * 10, Theta.thr = Theta.thr * 0.01,
                                  verbose = verbose, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                                  info = NULL, info.update = NULL, init.obj = init.obj, under.cv = FALSE)
      outTh.1se$Beta <- sweep(outTh.1se$Beta/init.obj$sdx, 2, init.obj$sdy, `*`)
      outTh.1se$mu <- as.numeric(init.obj$my - crossprod(outTh.1se$Beta, init.obj$mx))
      outTh.1se$lambda.Beta <- lamB.min
      outTh.1se$lambda.Theta <- lamTh.1se
      relax.net <- NULL
      if (fit.relax) {
        if (verbose > 0) { cat("    --- Relaxed network inference...\n\n") }
        relax.net <- relax.glasso(X = X, Y = Y, init.obj = init.obj, est = outTh.1se, eps = eps,
                                  Theta.thr = Theta.thr * 0.01, Theta.maxit = Theta.maxit * 10)
      }
      outTh.1se$relax.net <- relax.net
    } else {
      warning("\n`lambda.Theta.1se` = `lambda.Theta.min`, please provide a finer grid for `lambda.Theta` by increasing the number of values.\n")
    }
  }
  
  if (verbose > 0) { cat("========================= FINISHED ========================\n\n") }
  
  cv <- list(est.min = out.min, estB.1se = outB.1se, estTht.1se = outTh.1se, rho = init.obj$rho.vec, fold.index = ind,
             lambda.Beta.vec = lamB.vec, lambda.Theta.vec = lamTh.vec, cvm = err.cv, cvup = err.up, cvlo = err.low,
             penalize.diagonal = init.obj$penalize.diagonal, diag.penalty.factor = init.obj$diag.pf)
  class(cv) <- c("cv.missoNet", class(cv))
  return(cv)
}

