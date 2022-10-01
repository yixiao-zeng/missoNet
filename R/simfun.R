###################################
## Function for generating Beta*
###################################
genBeta <- function(p, q, s1, s2) {
  B <- matrix(rnorm(p * q, 0, 1), p, q)
  K <- matrix(0, p, q)
  R <- matrix(0, p, q)
  
  ## check the number of nonzero elements
  while ((sum(K * R) <= (p * q * s1 * s2)/1.5) | (sum(K * R) >= (p * q * s1 * s2) * 1.3)) {
    K <- matrix(rbinom(p * q, 1, s1), p, q)
    R <- matrix(0, p, q)
    tmp <- NULL
    while (length(tmp) == 0) {
      tmp <- which(rbinom(p, 1, s2) == 1)
    }
    R[tmp, ] <- 1
  }
  return(B * K * R)
}


###################################
## Function for generating Theta*
###################################
genTheta <- function(q) {
  Th <- matrix(0, q, q)
  while (!(min(eigen(Th)$value) > 1e-4)) {
    Th <- matrix(0, q, q)
    ## weak complete graph
    for (i in (round(q/4) + 1):(round(q/4) * 2)) {
      for (j in i:(round(q/4) * 2)) {
        if (i != j) {
          Th[i, j] <- Th[j, i] <- runif(1, min = -0.5, max = 0)
        }
      }
      Th[i, i] <- -sum(Th[i, (round(q/4) + 1):(round(q/4) * 2)])
    }
    ## strong complete graph
    for (i in (round(q/4) * 2 + 1):(round(q/4) * 3)) {
      for (j in i:(round(q/4) * 3)) {
        if (i != j) {
          Th[i, j] <- Th[j, i] <- runif(1, min = -1, max = -0.4)
        }
      }
      Th[i, i] <- -sum(Th[i, (round(q/4) * 2 + 1):(round(q/4) * 3)])
    }
    ## chain
    for (i in (round(q/4) * 3 + 1):q) {
      for (j in i:q) {
        if (i == (j - 1)) {
          Th[i, j] <- Th[j, i] <- runif(1, min = -1, max = -0.4)
        }
      }
      Th[i, i] <- -sum(Th[i, (round(q/4) * 3 + 1):q])
    }
    ## add small positive values to the diagonal elements to ensure PSD
    diag(Th) <- diag(Th) + 0.1
  }
  return(Th)
}


###################################
## Function for generating X
###################################
genX <- function(n, p, Sigma.X) {
  if (is.null(Sigma.X)) {
    s <- 0.7  # autocorrelation
    dist <- matrix(NA, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        dist[i, j] <- abs(i - j)
      }
    }
    Sigma.X <- s^dist
  } else {
    if (min(eigen(Sigma.X)$value) <= 1e-7) {
      stop("\nPlease supply a positive definite matrix for `Sigma.X`.\n")
    }
  }
  return(mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma.X, checkSymmetry = TRUE))
}


###################################
## Function for generating Y
###################################
genY <- function(X, BETAstar, Theta) {
  n <- nrow(X)
  q <- ncol(BETAstar)
  if (min(eigen(Theta)$value) <= 1e-7) {
    stop("\nPlease supply a positive definite matrix for `Theta`.\n")
  }
  Sigma <- solve(Theta)
  E <- mvtnorm::rmvnorm(n, rep(0, q), Sigma, checkSymmetry = TRUE)
  return(X %*% BETAstar + E)
}


###################################
## Function for generating Z
###################################
genZ <- function(X, BETAstar, Y, rho, type) {
  n <- nrow(Y)
  q <- ncol(Y)

  if (type == "MCAR") {
    missing <- matrix(1, n, q)
    for (j in 1:q) {
      missing[, j][rbinom(n, size = 1, prob = rho[j]) == 1] <- NA
    }
    Z <- Y * missing

  } else if (type == "MAR") {
    missing <- matrix(1, n, q)
    for (j in 1:q) {
      missing[, j][rbinom(n, size = 1, prob = 2 * rho[j] * 1/(1 + exp(-(X %*% BETAstar)[, j]))) == 1] <- NA
    }
    Z <- Y * missing

  } else if (type == "MNAR") {
    missing <- matrix(1, n, q)
    for (j in 1:q) {
      missing[, j][Y[, j] < quantile(Y[, j], rho[j])] <- NA
    }
    Z <- Y * missing
  }
  return(Z)
}

