############################################################################
# ADMM algorithm to create the positive definite estimate of the gram matrix.
############################################################################
maxproj.cov <- function(mat, epsilon = 1e-04, mu = 10, nitr.max = 1000, etol = 1e-04) {
  
  p <- nrow(mat)
  
  # Initialization
  R <- diag(mat)
  S <- matrix(0, p, p)
  L <- matrix(0, p, p)
  
  itr <- 0
  while (itr < nitr.max) {
    Rp <- R
    Sp <- S
    
    # Subproblem I: R step
    W <- mat + S + mu * L
    W.eigdec <- eigen(W, symmetric = TRUE)
    W.V <- W.eigdec$vectors
    W.D <- W.eigdec$values
    R <- W.V %*% diag(pmax(W.D, epsilon)) %*% t(W.V)
    
    # Subproblem II: S step
    M <- R - mat - mu * L
    S[lower.tri(S, diag = TRUE)] <- M[lower.tri(M, diag = TRUE)] - l1proj(v = M[lower.tri(M, diag = TRUE)], b = mu/2)
    for (i in 2:p) {
      for (j in 1:(i - 1)) {
        S[j, i] <- S[i, j]
      }
    }
    
    # L step: update the Lagrange parameter
    L <- L - (R - S - mat)/mu
    
    # Stopping Rule cat('check the stopping criterion:',max(abs(R-S-mat)),'\n')
    if ((max(abs(R - Rp)) < etol) && (max(abs(S - Sp)) < etol) && (max(abs(R - S - mat)) < etol)) {
      itr <- nitr.max
    } else {
      itr <- itr + 1
    }
    
    if (itr%%20 == 0) {
      mu <- mu/2
    }
  }
  return(R)
}


#######################################################################################################
# Efficient projection onto L1 ball of specified radius (i.e. b), used by the admm algo.
# Ref. Duchi et al. (2008). Efficient Projections onto the L1-Ball for Learning in High Dimensions, ICML.
#######################################################################################################
l1proj <- function(v, b) {
  
  stopifnot(b > 0)
  
  u <- sort(abs(v), decreasing = TRUE)
  sv <- cumsum(u)
  rho <- max(which(u > (sv - b)/1:length(u)))
  theta <- max(0, (sv[rho] - b)/rho)
  w <- sign(v) * pmax(abs(v) - theta, 0)
  
  return(w)
}


