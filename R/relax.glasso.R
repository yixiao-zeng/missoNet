relax.glasso <- function(X, Y, init.obj, est, eps, Theta.thr, Theta.maxit) {
  if (init.obj$penalize.diagonal) {
    thr <- 0.1
  } else {
    thr <- 0.01
  }
  
  Theta <- est$Theta
  q <- ncol(Theta)
  diag(Theta) <- 0
  if (max(abs(Theta)) != 0) {
    quant <- quantile(abs(unique(Theta[Theta!=0])), thr)
  } else {
    quant <- 0
  }
  
  zero <- NULL
  for (i in 1:(q - 1)) {
    for (j in (i + 1):q) {
      if (abs(Theta[i, j]) <= quant) {
        zero <- rbind(zero, c(i, j))
      }
    }
  }
  
  rho.mat <- matrix(1 - init.obj$rho.vec, q, 1) %*% matrix(1 - init.obj$rho.vec, 1, q)
  diag(rho.mat) <- 1 - init.obj$rho.vec
  n <- nrow(X)
  E <- Y - X %*% est$Beta
  E <- scale(E, scale = TRUE, center = TRUE)
  E[is.na(E)] <- 0
  residual.cov <- 1/n * (crossprod(E)/rho.mat)
  if (min(eigen(residual.cov)$value) < eps) {
    residual.cov <- maxproj.cov(mat = residual.cov, epsilon = eps)
  }
  # residual.cov <- getResidual(X = X, Y = Y, B = est$Beta, rho.mat = rho.mat, eps = eps)
  
  relax.graph <- suppressWarnings(glasso::glasso(s = residual.cov, rho = 0, thr = Theta.thr, maxit = Theta.maxit,
                                                 zero = zero, approx = FALSE, penalize.diagonal = TRUE, trace = FALSE))
  relax.graph <- (relax.graph$wi + t(relax.graph$wi))/2
  
  return(relax.graph)
}


