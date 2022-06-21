InitParams <- function(X, Y, rho, kfold, foldid, Theta.maxit, Theta.thr, eps, diag.pf, standardize, standardize.response) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  ## If n <= p (or q), the diagonal terms of Theta will be penalized
  penalize.diagonal <- (floor(n * (kfold - 1)/kfold) <= p | floor(n * (kfold - 1)/kfold) <= q)
  
  mx <- apply(X, 2, mean)
  my <- apply(Y, 2, mean, na.rm = TRUE)
  
  if (standardize) {
    sdx <- apply(X, 2, sd)
  } else {
    sdx <- rep(1, p)
  }
  
  if (standardize.response) {
    sdy <- apply(Y, 2, sd, na.rm = TRUE)
    Y <- scale(Y, center = FALSE, scale = sdy)
  } else {
    sdy <- rep(1, q)
  }
  
  B.init <- matrix(0, nrow = p, ncol = q)
  lam.list <- rep(0, q)
  for (j in 1:q) {
    cc <- complete.cases(Y[, j])
    cv <- glmnet::cv.glmnet(x = X[cc, ], y = Y[cc, j], intercept = TRUE, standardize = standardize, 
                            family = "gaussian", foldid = foldid[cc])
    B.init[, j] <- coef(cv, s = cv$lambda.min)[2:(p + 1)]
    lam.list[j] <- max(cv$lambda)
  }
  lamB.max <- max(lam.list) * 4
  
  
  if (is.null(rho)) {
    rho.vec <- apply(Y, 2, function(x) {
      sum(is.na(x))/n
    })
  } else if (length(rho) == 1) {
    rho.vec <- rep(rho, q)
  } else if (length(rho) == q) {
    rho.vec <- rho
  } else {
    stop("rho must be a scalar or a vector of length = q.")
  }
  
  rho.mat <- matrix(1 - rho.vec, q, 1) %*% matrix(1 - rho.vec, 1, q)
  diag(rho.mat) <- 1 - rho.vec  # qxq
  
  residual.cov.init <- getResidual(X = X, Y = Y, B = B.init, rho.mat = rho.mat, eps = eps)
  glasso.obj <- glasso::glassopath(s = residual.cov.init, rholist = NULL, thr = Theta.thr, maxit = Theta.maxit, 
                                   approx = FALSE, penalize.diagonal = penalize.diagonal, trace = 0)
  if (penalize.diagonal) {
    if (is.null(diag.pf)) {
      diag.wi <- diag(diag(glasso.obj$wi[,,1]))
      offdiag.wi <- glasso.obj$wi[,,1] - diag.wi
      if(max(abs(offdiag.wi)) != 0) {
        diag.pf <- max(abs(diag.wi))/max(abs(offdiag.wi))
      }else{
        diag.pf <- 4
        warning("Unable to compute diag.penalty.factor, set to the default value = 4.
Supply different values to tune results.")
      }
    }
  } else {
    diag.pf <- NULL
  }
  lamTh.max <- max(glasso.obj$rholist) * 2
  
  return(list(B.init = B.init, residual.cov.init = residual.cov.init, rho.vec = rho.vec, sdx = sdx, sdy = sdy, mx = mx, my = my,
              penalize.diagonal = penalize.diagonal, diag.pf = diag.pf, lamB.max = lamB.max, lamTh.max = lamTh.max))
}


