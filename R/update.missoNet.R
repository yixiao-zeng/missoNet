update.missoNet <- function(X, Y, lamTh, lamB,
                            Beta.maxit, Beta.thr, Theta.maxit, Theta.thr,
                            verbose, eps, eta, diag.pf,
                            info, info.update, under.cv, init.obj=NULL, B.init=NULL) {
  if (is.null(info)) {
    n <- nrow(X)
    p <- ncol(X)
    q <- ncol(Y)
    
    X <- scale(X, center = init.obj$mx, scale = init.obj$sdx)
    Y <- scale(Y, center = init.obj$my, scale = init.obj$sdy)
    Z <- Y
    Z[is.na(Z)] <- 0
    
    rho.mat.1 <- t(matrix(rep(1 - init.obj$rho.vec, p), q, p))  # pxq
    rho.mat.2 <- matrix(1 - init.obj$rho.vec, q, 1) %*% matrix(1 - init.obj$rho.vec, 1, q)
    diag(rho.mat.2) <- 1 - init.obj$rho.vec  # qxq
    
    info$n <- n
    info$q <- q
    info$penalize.diagonal <- init.obj$penalize.diagonal
    info$xtx <- crossprod(X)
    info$til.xty <- crossprod(X, Z)/rho.mat.1
    til.ytx <- t(info$til.xty)
    til.yty <- crossprod(Z)/rho.mat.2
    if (min(eigen(til.yty)$value) < eps) {
      til.yty <- maxproj.cov(mat = til.yty, epsilon = eps)
    }
    
    #####################################################
    # Pre-updating several times due to a cold start
    #####################################################
    if (verbose == 2) {
      cat("  ---------------- Warming-up -----------------\n")
      cat("\titer\t|\t| lik(t + 1) - lik(t) |\n")
    }
    Beta.thr.rescale <- Beta.thr * sum(abs(B.init))
    E <- Y - X %*% B.init
    residual.cov <- getResidual(E = E, n = n, rho.mat = rho.mat.2, eps = eps)
    
    if (info$penalize.diagonal) {
      lamTh.mat <- lamTh * (1 - diag(info$q)) + lamTh * diag.pf * diag(info$q)
    } else {
      lamTh.mat <- lamTh * (1 - diag(info$q))
    }
    lamB.mat <- matrix(lamB, nrow = p, ncol = q)
    
    lik.new <- sum(diag(1/n * (til.yty - til.ytx %*% B.init - crossprod(B.init, info$til.xty)
                               + crossprod(B.init, info$xtx) %*% B.init) %*% diag(1, q)))
    - determinant(diag(1, q), logarithm = TRUE)$mod[1] + sum(abs(lamTh.mat * diag(1, q))) + sum(abs(lamB.mat * B.init))
    lik.thr <- lik.new * 1e-07
    lik.old <- lik.new + lik.thr + 1
    
    s <- 0
    while(s < 50) {
      if(abs(lik.new - lik.old) < lik.thr) {
        if (verbose == 2) {
          cat("  ---------------------------------------------\n")
        }
        break
      } else {
        lik.old <- lik.new
        Theta.out <- glasso(s = residual.cov, rho = lamTh.mat, thr = Theta.thr, maxit = Theta.maxit,
                            approx = FALSE, penalize.diagonal = info$penalize.diagonal, trace = FALSE)
        Theta <- (Theta.out$wi + t(Theta.out$wi))/2

        B.out <- updateBeta(Theta = Theta, B0 = B.init, n = info$n, xtx = info$xtx, xty = info$til.xty,
                            lamB = lamB, eta = eta, tolin = Beta.thr.rescale, maxitrin = Beta.maxit)
        
        B.init <- B.out$Bhat
        Beta.thr.rescale <- Beta.thr * sum(abs(B.init))
        E <- Y - X %*% B.init
        residual.cov <- getResidual(E = E, n = n, rho.mat = rho.mat.2, eps = eps)
        
        lik.new <- sum(diag(1/n * (til.yty - til.ytx %*% B.init - crossprod(B.init, info$til.xty)
                                   + crossprod(B.init, info$xtx) %*% B.init) %*% Theta))
        - determinant(Theta, logarithm = TRUE)$mod[1] + sum(abs(lamTh.mat * Theta)) + sum(abs(lamB.mat * B.init))
        
        if (verbose == 2) {
          cat(sprintf("\t%d\t|\t%f\n", (s + 1), abs(lik.new - lik.old)))
        }
        s <- s + 1
      }
    }
    #####################################################
    # Pre-updating ends
    #####################################################
    info.update$B.init <- B.init
    Beta.thr <- Beta.thr.rescale
    info.update$residual.cov <- residual.cov
  }
  
  ################################################################################
  # Updating Theta and Beta
  ################################################################################
  if (info$penalize.diagonal) {
    lamTh.mat <- lamTh * (1 - diag(info$q)) + lamTh * diag.pf * diag(info$q)
    Theta.out <- glasso(s = info.update$residual.cov, rho = lamTh.mat, thr = Theta.thr, maxit = Theta.maxit,
                        approx = FALSE, penalize.diagonal = TRUE, trace = FALSE)
  } else {
    Theta.out <- glasso(s = info.update$residual.cov, rho = lamTh, thr = Theta.thr, maxit = Theta.maxit,
                        approx = FALSE, penalize.diagonal = FALSE, trace = FALSE)
  }
  Theta <- (Theta.out$wi + t(Theta.out$wi))/2
  
  B.out <- updateBeta(Theta = Theta, B0 = info.update$B.init, n = info$n, xtx = info$xtx, xty = info$til.xty,
                      lamB = lamB, eta = eta, tolin = Beta.thr, maxitrin = Beta.maxit)
  
  if (verbose == 2) {
    cat("  `lambda.Theta`:", lamTh, "   `lambda.Beta`:", lamB, "\n")
    cat("  # iters for updating `Theta`: ", Theta.out$niter, "\n")
    cat("  # iters for updating `Beta`: ", B.out$it.final, "\n\n")
  }
  
  if(under.cv) {
    return(B.out$Bhat)
  } else {
    return(list(Beta = B.out$Bhat, Theta = Theta))
  }
}

