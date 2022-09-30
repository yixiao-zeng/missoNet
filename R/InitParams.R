InitParams <- function(X, Y, rho, kfold, foldid, Theta.maxit, Theta.thr, eps,
                       penalize.diagonal, diag.pf, standardize, standardize.response) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  ## If n.tr <= max(p, q), the diagonal of Theta should be penalized
  if (is.null(penalize.diagonal)) {
    penalize.diagonal <- (floor(n * (kfold - 1)/kfold) <= max(p, q))
  } else {
    if ((floor(n * (kfold - 1)/kfold) <= max(p, q)) & (penalize.diagonal == FALSE)) {
      warning("\nInsufficient sample size. It is recommended to set `penalize.diagonal = TRUE`.\n")
    } else if ((floor(n * (kfold - 1)/kfold) > max(p, q)) & (penalize.diagonal == TRUE)) {
      warning("\nSufficient sample size. It is recommended to set `penalize.diagonal = FALSE`.\n")
    }
  }
  
  mx <- apply(X, 2, mean)
  my <- apply(Y, 2, mean, na.rm = TRUE)
  
  if (standardize) {
    sdx <- apply(X, 2, sd)
  } else { sdx <- rep(1, p) }
  
  if (standardize.response) {
    sdy <- apply(Y, 2, sd, na.rm = TRUE)
    Y <- scale(Y, center = FALSE, scale = sdy)
  } else { sdy <- rep(1, q) }
  
  B.init <- matrix(0, nrow = p, ncol = q)
  lam.list <- rep(0, q)
  for (j in 1:q) {
    cc <- complete.cases(Y[, j])
    cv <- glmnet::cv.glmnet(x = X[cc, ], y = Y[cc, j], intercept = TRUE, standardize = standardize, 
                            family = "gaussian", foldid = foldid[cc])
    B.init[, j] <- coef(cv, s = cv$lambda.min)[2:(p + 1)]
    lam.list[j] <- max(cv$lambda)
  }
  lamB.max <- max(lam.list) * 10
  
  
  if (is.null(rho)) {
    rho.vec <- apply(Y, 2, function(x) {
      sum(is.na(x))/n
    })
  } else if (length(rho) == 1) {
    rho.vec <- rep(rho, q)
  } else if (length(rho) == q) {
    rho.vec <- rho
  } else {
    stop("\n`rho` must be a scalar or a vector of length = q.\n")
  }
  
  rho.mat <- matrix(1 - rho.vec, q, 1) %*% matrix(1 - rho.vec, 1, q)
  diag(rho.mat) <- 1 - rho.vec  # qxq
  
  E <- Y - X %*% B.init
  residual.cov.init <- getResidual(E = E, n = n, rho.mat = rho.mat, eps = eps)
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
        warning("\nFail to compute the `diag.penalty.factor`, has been set to the default value = 4.\n
Please try to reduce the dimensionality by filtering variables prior to fitting the model, 
or try different values (usually larger) for the argument `diag.penalty.factor`.\n")
      }
    }
  } else {
    diag.pf <- NULL
  }
  lamTh.max <- max(glasso.obj$rholist) * 10
  
  return(list(B.init = B.init, residual.cov.init = residual.cov.init, rho.vec = rho.vec, sdx = sdx, sdy = sdy, mx = mx, my = my,
              penalize.diagonal = penalize.diagonal, diag.pf = diag.pf, lamB.max = lamB.max, lamTh.max = lamTh.max))
}

