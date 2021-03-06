#' Cross-validation for missoNet
#'
#' @param X Numeric predictor matrix (n by p): columns correspond to predictor variables and rows correspond to samples. Missing values are not allowed. Do not include a column of ones.
#' @param Y Numeric response matrix (n by q): columns correspond to response variables and rows correspond to samples. Missing values should be coded as \code{NA} or \code{NaN}.
#' @param kfold Number of folds -- default is 5.
#' @param rho (Optional) A scalar or a numeric vector of length q: a user supplied missing probability for response variables. Default is \code{rho = NULL} and the program will compute the empirical missing rates for columns of \code{Y} and use them as the missing probability.
#' @param lambda.Beta (Optional) Numeric vector: a user supplied regularization parameter sequence for \code{Beta} from which the CV procedure searches; default is \code{lambda.Beta = NULL} and the program computes its own lambda sequence based on \code{n.lamBeta} and \code{lamBeta.min.ratio}, supplying a vector overrides this. WARNING: avoid supplying a single value, missoNet relies on warm starts for acceleration and accuracy. The sequence supplied will be automatically arranged in descending order internally. 
#' @param lambda.Theta (Optional) Numeric vector: a user supplied regularization parameter sequence for \code{Theta} from which the CV procedure searches; default is \code{lambda.Theta = NULL} and the program computes its own lambda sequence based on \code{n.lamTheta} and \code{lamTheta.min.ratio}, supplying a vector overrides this. WARNING: avoid supplying a single value, missoNet relies on warm starts for acceleration and accuracy. The sequence supplied will be automatically arranged in descending order internally. 
#' @param lamBeta.min.ratio The smallest value for \code{lambda.Beta} will be data-derived \code{max(lambda.Beta)} multiplied by \code{lamBeta.min.ratio}. The default depends on the sample size n relative to the number of predictors p. If n > p, default is 0.01; if n <= p, default is 0.1. A very small value of \code{lamBeta.min.ratio} will lead to long computation time and may cause a saturated fit in the n < p case. Only needed when \code{lambda.Beta = NULL}.
#' @param lamTheta.min.ratio The smallest value for \code{lambda.Theta} will be data-derived \code{max(lambda.Theta)} multiplied by \code{lamTheta.min.ratio}; default is 0.01. A very small value of \code{lamTheta.min.ratio} will lead to long computation time. Only needed when \code{lambda.Theta = NULL}.
#' @param lamBeta.scale.factor A multiplication factor for scaling the entire \code{lambda.Beta} sequence; defaults is 1. A typical usage scenario is when the optimal lambda approaches the boundary of the search area. Only needed when \code{lambda.Beta = NULL}.
#' @param lamTheta.scale.factor A multiplication factor for scaling the entire \code{lambda.Theta} sequence; defaults is 1. A typical usage scenario is when the optimal lambda approaches the boundary of the search area. Only needed when \code{lambda.Theta = NULL}.
#' @param n.lamBeta The number of \code{lambda.Beta} values. By default, the program estimates a number based on \code{lamBeta.min.ratio}. An excessively large number results in a significant increase in computation time. Only needed when \code{lambda.Beta = NULL}.
#' @param n.lamTheta The number of \code{lambda.Theta} values. By default, the program estimates a number based on \code{lamTheta.min.ratio}. An excessively large number results in a significant increase in computation time. Only needed when \code{lambda.Theta = NULL}.
#' @param Beta.maxit The maximum number of internal iterations allowed for updating \code{Beta}. Default is \code{Beta.maxit = 1000}.
#' @param Beta.thr The convergence threshold for updating \code{Beta}; default is \code{Beta.thr = 1e-4}. Iterations stop when absolute parameter change is less than \code{Beta.thr * sum(abs(Beta))}.
#' @param eta Backtracking line search shrinkage factor, default is \code{eta = 0.8}. You may want to choose a more appropriate \code{eta} for a faster \code{Beta} convergence based on your dataset. Note that \code{eta} must be greater than 0 and smaller than 1.
#' @param Theta.maxit The maximum number of internal iterations allowed for updating \code{Theta}. Default is \code{Theta.maxit = 1000}.
#' @param Theta.thr The convergence threshold for updating \code{Theta}; default is \code{Theta.thr = 1e-4}. Iterations stop when average absolute parameter change is less than \code{Theta.thr * ave(abs(offdiag(Sigma)))}.
#' @param eps A numeric tolerance level for L1 projection; default is \code{eps = 1e-8}. If any of the eigenvalues is less than the given tolerance, the unbiased estimate of covariance is projected onto L1 ball to have \code{min(eigen(Sigma)$value) == eps}.
#' @param diag.penalty.factor Numeric: a separate penalty factor for the diagonal entries of \code{Theta}. This is a number that multiplies \code{lambda.Theta} to allow differential shrinkage. Default is \code{NULL} and the program computes it based on an initial estimate of \code{Theta}. Can be 0 for no shrinkage. Only needed when n <= p.
#' @param standardize Logical: should the columns of \code{X} be standardized so each has unit length and zero average; default is \code{TRUE}. The parameter estimates will be returned on the original scale. If \code{X} has been standardized prior to fitting the model, you might not wish to standardize.
#' @param standardize.response Logical: should the columns of \code{Y} be standardized so each has unit length and zero average; default is \code{TRUE}. The parameter estimates will be returned on the original scale. If \code{Y} has been standardized prior to fitting the model, you might not wish to standardize.
#' @param fit.1se Logical: default is \code{FALSE}. Should the model be refitted with the largest \code{lambda.Beta} according to the one-standard-error rule?
#' @param fit.relax Logical: default is \code{FALSE}. If \code{TRUE}, the program will re-estimate a relaxed graph (\code{Theta}) without penalization, which could be useful for network analysis. WARNING: there may be convergence issues if the residual covariance matrix is not of full rank.
#' @param permute Logical: should the subject indices for CV be permuted? Default is \code{FALSE}.
#' @param with.seed A random seed for permutation.
#' @param parallel Logical: default is \code{FALSE}. If \code{TRUE}, use parallel cluster to fit each fold.
#' @param cpus Number of cores for parallelization. Only needed when \code{parallel = TRUE}.
#' @param verbose Value of 0, 1 or 2. 0 -- silent; 1 -- limited tracing; 2 -- detailed tracing. Limited tracing if \code{parallel = TRUE}.
#'
#' @return This function returns a cv.missoNet object containing a \code{list} of estimates and working parameters:
#' \itemize{
#' \item \code{est.min}: estimates using the pair of regularization parameters that gives the minimum cross-validation error.
#' \item \code{est.1se}: estimates using the pair of regularization parameters such that error is within one-standard-error of the minimum, if \code{fit.1se = TRUE}.
#' \item \code{lambda.Beta.vec}: a long vector of length \code{n.lamBeta * n.lamTheta} from which the CV procedure searches.
#' \item \code{lambda.Theta.vec}: a long vector of length \code{n.lamBeta * n.lamTheta} from which the CV procedure searches.
#' \item \code{relax.graph}: an estimate of network without penalization if \code{fit.relax = TRUE}.
#' \item \code{rho}: a vector of length q: working missing probability.
#' \item \code{fold.index}: subject indices for cross-validation split.
#' \item \code{cvm}: mean cross-validated error that has a one to one correspondence with \code{lambda.Beta.vec} and \code{lambda.Theta.vec}.
#' \item \code{cvup}: upper cross-validated error.
#' \item \code{cvlo}: lower cross-validated error.
#' \item \code{penalize.diagonal}: logical: whether the diagonal of \code{Theta} was penalized.
#' \item \code{diag.penalty.factor}: penalty factor for the diagonal when \code{penalize.diagonal} returned as \code{TRUE}.
#' }
#' @export
#' 
#' @examples 

cv.missoNet <- function(X, Y, kfold = 5, rho = NULL, lambda.Beta = NULL, lambda.Theta = NULL,
                        lamBeta.min.ratio = 0.01, lamTheta.min.ratio = 0.01, lamBeta.scale.factor = 1, lamTheta.scale.factor = 1,
                        n.lamBeta = round(-log10(lamBeta.min.ratio) * 15), n.lamTheta = round(-log10(lamTheta.min.ratio) * 10),
                        Beta.maxit = 1000, Beta.thr = 1e-04, eta = 0.8, Theta.maxit = 1000, Theta.thr = 1e-04, eps = 1e-08, 
                        diag.penalty.factor = NULL, standardize = TRUE, standardize.response = TRUE, fit.1se = FALSE, fit.relax = FALSE,
                        permute = FALSE, with.seed = NULL, parallel = FALSE, cpus, verbose = 1) {
  if (verbose > 0) {
    cat("\nInitializing necessary parameters...\n\n")
  }
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  if (permute) {
    set.seed(with.seed)
    ind <- sample(n, replace = FALSE)
  } else {
    ind <- 1:n
  }
  foldid <- unlist(lapply(1:kfold, function(x) {
    rep(x, length((1 + floor((x - 1) * n/kfold)):floor(x * n/kfold)))
  }))
  names(ind) <- paste0("fold.", foldid)
  
  init.obj <- InitParams(X = X[ind, ], Y = Y[ind, ], rho = rho, kfold = kfold, foldid = foldid,
                         Theta.maxit = Theta.maxit, Theta.thr = Theta.thr, eps = eps, diag.pf = diag.penalty.factor,
                         standardize = standardize, standardize.response = standardize.response)
  
  lambda.obj <- InitLambda(lamB = lambda.Beta, lamTh = lambda.Theta, 
                           init.obj = init.obj, n.lamB = n.lamBeta, n.lamTh = n.lamTheta,
                           lamB.min.ratio = lamBeta.min.ratio, lamTh.min.ratio = lamTheta.min.ratio,
                           lamB.scale.factor = lamBeta.scale.factor, lamTh.scale.factor = lamTheta.scale.factor)
  lamTh.vec <- lambda.obj$lamTh.vec
  lamB.vec <- lambda.obj$lamB.vec
  
  ################################################################################
  # Cross-validation
  ################################################################################
  if (!parallel) {
    err <- matrix(NA, kfold, length(lamTh.vec))
    for (k in 1:kfold) {
      if (verbose > 0) {
        cat("Staring fold: ", k, "\n")
      }
      
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
      } else {
        rho.vec <- init.obj$rho.vec
      }
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
      }
    }
    
  } else {
    if (verbose > 0) {
      snowfall::sfInit(parallel = TRUE, cpus = cpus)
    } else {
      suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = cpus))
    }
    
    par.out <- snowfall::sfClusterApplyLB(seq(1, kfold, by = 1), function(k) {
      parWrapper(k = k, X = X, Y = Y, init.obj = init.obj, rho = rho, ind = ind, kfold = kfold, lamTh.vec = lamTh.vec, lamB.vec = lamB.vec,
                 Beta.maxit = Beta.maxit, Beta.thr = Beta.thr, Theta.maxit = Theta.maxit, Theta.thr = Theta.thr, eps = eps, eta = eta)
    })
    
    if (verbose > 0) {
      snowfall::sfStop()
    } else {
      suppressMessages(snowfall::sfStop())
    }
    
    err <- do.call("rbind", par.out)
  }
  
  err.cv <- colSums(err)/kfold
  err.sd <- apply(err, 2, sd)/sqrt(kfold)
  err.up <- err.cv + err.sd
  err.low <- err.cv - err.sd
  
  cv.min <- which.min(err.cv)
  lamTh.min <- lamTh.vec[cv.min]
  lamB.min <- lamB.vec[cv.min]
  if (lamTh.min == max(lamTh.vec)) {
    warning("Optimal lambda.Theta is at the boundary, provide a larger value to one of the following:
    1. lamTheta.scale.factor;
    2. lambda.Theta (if user specified).\n")
  } else if (lamTh.min == min(lamTh.vec)) {
    warning("Optimal lambda.Theta is at the boundary, provide a smaller value to one of the following:
    1. lamTheta.scale.factor;
    2. lamTheta.min.ratio;
    3. lambda.Theta (if user specified).\n")
  }
  if (lamB.min == max(lamB.vec)) {
    warning("Optimal lambda.Beta is at the boundary, provide a larger value to one of the following:
    1. lamBeta.scale.factor;
    2. lambda.Beta (if user specified).\n")
  } else if (lamB.min == min(lamB.vec)) {
    warning("Optimal lambda.Beta is at the boundary, provide a smaller value to one of the following:
    1. lamBeta.scale.factor;
    2. lamBeta.min.ratio;
    3. lambda.Beta (if user specified).\n")
  }
  
  if (verbose > 0) {
    cat("Cross validation completed.\n\n")
    cat("Fittig with the lambda pair that gives the minimum CV error...\n\n")
  }
  out.min <- update.missoNet(X = X, Y = Y, lamTh = lamTh.min, lamB = lamB.min,
                             Beta.maxit = Beta.maxit * 10, Beta.thr = Beta.thr * 0.1,
                             Theta.maxit = Theta.maxit * 10, Theta.thr = Theta.thr * 0.1,
                             verbose = verbose, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                             info = NULL, info.update = NULL, init.obj = init.obj, under.cv = FALSE)
  out.min$lambda.Beta <- lamB.min
  out.min$lambda.Theta <- lamTh.min
  out.min$Beta <- sweep(out.min$Beta/init.obj$sdx, 2, init.obj$sdy, `*`)    ## convert back to the original scale
  out.min$mu <- as.numeric(init.obj$my - crossprod(out.min$Beta, init.obj$mx))
  if (verbose > 0) {
    cat("Done.\n\n")
  }
  
  out.1se <- NULL
  if (fit.1se) {
    if (verbose > 0) {
      cat("Fiting with the largest lambda.Beta such that error is within 1se of the minimum....\n\n")
    }
    new.lamB.vec <- lamB.vec[lamTh.vec == lamTh.min]
    new.err.cv <- err.cv[lamTh.vec == lamTh.min]
    tmp <- new.lamB.vec >= lamB.min
    new.lamB.vec <- new.lamB.vec[tmp]
    new.err.cv <- new.err.cv[tmp]
    lamB.1se <- new.lamB.vec[which.min(abs(new.err.cv - err.up[cv.min]))]
    
    out.1se <- update.missoNet(X = X, Y = Y, lamTh = lamTh.min, lamB = lamB.1se,
                               Beta.maxit = Beta.maxit * 10, Beta.thr = Beta.thr * 0.1,
                               Theta.maxit = Theta.maxit * 10, Theta.thr = Theta.thr * 0.1,
                               verbose = verbose, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                               info = NULL, info.update = NULL, init.obj = init.obj, under.cv = FALSE)
    out.1se$lambda.Beta <- lamB.1se
    out.1se$lambda.Theta <- lamTh.min
    out.1se$Beta <- sweep(out.1se$Beta/init.obj$sdx, 2, init.obj$sdy, `*`)
    out.1se$mu <- as.numeric(init.obj$my - crossprod(out.1se$Beta, init.obj$mx))
    if (verbose > 0) {
      cat("Done.\n\n")
    }
  }
  
  relax.graph <- NULL
  if (fit.relax) {
    relax.graph <- relax.glasso(X = X, Y = Y, init.obj = init.obj, est = out.min, eps = eps,
                                Theta.thr = Theta.thr * 0.1, Theta.maxit = Theta.maxit * 10)
  }
  
  cv <- list(est.min = out.min, est.1se = out.1se, lambda.Beta.vec = lamB.vec, lambda.Theta.vec = lamTh.vec,
             relax.graph = relax.graph, rho = init.obj$rho.vec, fold.index = ind, cvm = err.cv, cvup = err.up, cvlo = err.low,
             penalize.diagonal = init.obj$penalize.diagonal, diag.penalty.factor = init.obj$diag.pf)
  class(cv) <- c("cv.missoNet", class(cv))
  return(cv)
}


