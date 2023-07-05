InitParams <- function(X, Y, rho, under.cv, lamB.vec,
                       eps, penalize.diag, diag.pf,
                       standardize, standardize.response) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  ## Standardization
  mx <- apply(X, 2, mean)
  if (standardize) {
    sdx <- apply(X, 2, sd)
  } else { sdx <- rep(1, p) }
  X <- scale(X, center = mx, scale = sdx)
  
  my <- apply(Y, 2, mean, na.rm = TRUE)
  if (standardize.response) {
    sdy <- apply(Y, 2, sd, na.rm = TRUE)
  } else { sdy <- rep(1, q) }
  Y <- scale(Y, center = my, scale = sdy)
  Z <- Y
  Z[is.na(Z)] <- 0
  
  ## Initialization of un-biased surrogates
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
  
  rho.mat.1 <- t(matrix(rep(1 - rho.vec, p), q, p))  # pxq
  rho.mat.2 <- matrix(1 - rho.vec, q, 1) %*% matrix(1 - rho.vec, 1, q)
  diag(rho.mat.2) <- 1 - rho.vec  # qxq
  
  xtx <- crossprod(X)
  til.xty <- crossprod(X, Z)/rho.mat.1
  
  ## B.init and lambda.Beta.max
  if (is.null(lamB.vec)) {
    lamB.max <- 2/n * max(abs(til.xty)) * 5
  } else {
    lamB.max <- max(lamB.vec)
  }
  
  if (under.cv) {
    B.init <- updateBeta(Theta = diag(1, q), B0 = matrix(0, nrow = p, ncol = q), n = n, xtx = xtx, xty = til.xty,
                         lamB = lamB.max, eta = 0.8, tolin = 1e-08, maxitrin = 1e3)$Bhat
  } else {
    B.init <- lapply(lamB.vec, function(b) {
      updateBeta(Theta = diag(1, q), B0 = matrix(0, nrow = p, ncol = q), n = n, xtx = xtx, xty = til.xty,
                 lamB = b, eta = 0.8, tolin = 1e-08, maxitrin = 1e3)$Bhat })
  }
  
  ## lambda.Theta.max
  residual.cov <- 1/n * (crossprod(Z)/rho.mat.2)
  if (min(eigen(residual.cov)$value) < eps) { residual.cov <- maxproj.cov(mat = residual.cov, epsilon = eps) }
  lamTh.max <- max(abs(residual.cov)) * 10
  if (penalize.diag) {
    if (is.null(diag.pf)) {
      glasso.obj <- suppressWarnings(glasso::glasso(s = residual.cov, rho = 0, approx = FALSE, penalize.diagonal = TRUE, trace = FALSE))
      diag.w <- diag(glasso.obj$w)
      offdiag.w <- glasso.obj$w[lower.tri(glasso.obj$w) | upper.tri(glasso.obj$w)]
      if (max(abs(offdiag.w)) != 0) {
        diag.pf <- max(max(abs(diag.w))/max(abs(offdiag.w)), 1)
      } else {
        diag.pf <- 1
        warning("\nFail to compute the `diag.penalty.factor`, has been set to the default value = 1.\n
Please try to reduce the dimensionality by filtering variables prior to fitting the model,
or try (usually) larger values for the argument `diag.penalty.factor`.\n")
      }
    }
  } else {
    diag.pf <- NULL
  }

  return(list(B.init = B.init, rho.vec = rho.vec, sdx = sdx, sdy = sdy, mx = mx, my = my, 
              lamB.max = lamB.max, lamTh.max = lamTh.max, diag.pf = diag.pf))
}

