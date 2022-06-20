getResidual <- function(X, Y, B, rho.mat, eps) {
  n <- nrow(X)
  E <- Y - X %*% B
  E <- scale(E, scale = FALSE, center = TRUE)
  E[is.na(E)] <- 0
  residual.cov <- 1/n * (crossprod(E)/rho.mat)
  if (min(eigen(residual.cov)$value) < eps) {
    residual.cov <- maxproj.cov(mat = residual.cov, epsilon = eps)
  }
  
  return(residual.cov)
}

