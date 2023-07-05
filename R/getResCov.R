getResCov <- function(E, n, rho.mat, eps) {
  E <- scale(E, scale = FALSE, center = TRUE)
  E[is.na(E)] <- 0
  residual.cov <- 1/n * (crossprod(E)/rho.mat)
  if (min(eigen(residual.cov)$value) < eps) {
    residual.cov <- maxproj.cov(mat = residual.cov, epsilon = eps)
  }
  return(residual.cov)
}

# getResCov <- function(yty, ytx, xty, xtx, B, n, eps) {
#   residual.cov <- 1/n * (yty - ytx %*% B - crossprod(B, xty) + crossprod(B, xtx) %*% B)
#   if (min(eigen(residual.cov)$value) < eps) {
#     residual.cov <- maxproj.cov(mat = residual.cov, epsilon = eps)
#   }
#   return(residual.cov)
# }
