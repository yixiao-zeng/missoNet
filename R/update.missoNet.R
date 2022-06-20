update.missoNet <- function(X, Y, lamTh, lamB,
                            Beta.maxit, Beta.thr, Theta.maxit, Theta.thr,
                            verbose, eps, diag.pf,
                            info, init.obj) {
  if (is.null(info)) {
    n <- nrow(X)
    p <- ncol(X)
    q <- ncol(Y)
    
    X <- scale(X, center = init.obj$mx, scale = init.obj$sdx)
    Y <- scale(Y, center = init.obj$my, scale = init.obj$sdy)
    Y[is.na(Y)] <- 0
    
    rho.mat.1 <- t(matrix(rep(1 - init.obj$rho.vec, p), q, p))  # pxq
    rho.mat.2 <- matrix(1 - init.obj$rho.vec, q, 1) %*% matrix(1 - init.obj$rho.vec, 1, q)
    diag(rho.mat.2) <- 1 - init.obj$rho.vec  # qxq
    
    info <- NULL
    info$n <- n
    info$q <- q
    info$penalize.diagonal <- init.obj$penalize.diagonal
    info$xtx <- crossprod(X)
    info$til.xty <- crossprod(X, Y)/rho.mat.1
    til.ytx <- t(info$til.xty)
    til.yty <- crossprod(Y)/rho.mat.2
    if (min(eigen(til.yty)$value) < eps) {
      til.yty <- maxproj.cov(mat = til.yty, epsilon = eps)
    }
    
    #####################################################
    # Pre-updating several times due to a non-warming start
    #####################################################
    B.init <- init.obj$B.init * init.obj$sdx
    Beta.thr.rescale <- Beta.thr * sum(abs(B.init))
    residual.cov <- init.obj$residual.cov.init
    # w.init <- matrix(0, nrow = q, ncol = q)
    # wi.init <- matrix(0, nrow = q, ncol = q)
    
    if (info$penalize.diagonal) {
      lamTh.mat <- lamTh * (1 - diag(info$q)) + lamTh * diag.pf * diag(info$q)
    } else {
      lamTh.mat <- lamTh * (1 - diag(info$q))
    }
    lamB.mat <- matrix(lamB, nrow = p, ncol = q)
    
    lik.new <- sum(diag(1/n * (til.yty - til.ytx %*% B.init - crossprod(B.init, info$til.xty)
                               + crossprod(B.init, info$xtx) %*% B.init) %*% diag(1, q)))
    - determinant(diag(1, q), logarithm = TRUE)$mod[1] + sum(abs(lamTh.mat * diag(1, q))) + sum(abs(lamB.mat * B.init))
    lik.thr <- lik.new * 1e-06
    lik.old <- lik.new + lik.thr + 1
    
    s <- 0
    while(s < 50) {
      if(abs(lik.new - lik.old) < lik.thr) {
        if (verbose == 2) {
          cat("  Pre-updating completed.\n\n")
        }
        break
      } else {
        lik.old <- lik.new
        Theta.out <- glasso::glasso(s = residual.cov, rho = lamTh.mat, thr = Theta.thr, maxit = Theta.maxit,
                                    approx = FALSE, penalize.diagonal = info$penalize.diagonal, trace = FALSE)
        Theta <- (Theta.out$wi + t(Theta.out$wi))/2

        B.out <- updateBeta(Theta = Theta, B0 = B.init, lamB = lamB, eta = 0.8,
                            tolin = Beta.thr.rescale, maxitrin = Beta.maxit, info = info)
        
        B.init <- B.out$Bhat
        Beta.thr.rescale <- Beta.thr * sum(abs(B.init))
        residual.cov <- getResidual(X = X, Y = Y, B = B.init, rho.mat = rho.mat.2, eps = eps)
        # w.init <-  Theta.out$w
        # wi.init <-  Theta
        
        lik.new <- sum(diag(1/n * (til.yty - til.ytx %*% B.init - crossprod(B.init, info$til.xty)
                                   + crossprod(B.init, info$xtx) %*% B.init) %*% Theta))
        - determinant(Theta, logarithm = TRUE)$mod[1] + sum(abs(lamTh.mat * Theta)) + sum(abs(lamB.mat * B.init))
        
        if (verbose == 2) {
          cat("  Pre-updating iteration ", (s + 1), " -- likelihood change: ", abs(lik.new - lik.old), "\n")
        }
        s <- s + 1
      }
    }
    #####################################################
    # Pre-updating ends
    #####################################################
    
    info$B.init <- B.init
    Beta.thr <- Beta.thr.rescale
    info$residual.cov <- residual.cov
    # info$w.init <- w.init
    # info$wi.init <- wi.init
  }
  
  ################################################################################
  # Updating Theta and Beta
  ################################################################################
  if (info$penalize.diagonal) {
    lamTh.mat <- lamTh * (1 - diag(info$q)) + lamTh * diag.pf * diag(info$q)
    Theta.out <- glasso::glasso(s = info$residual.cov, rho = lamTh.mat, thr = Theta.thr, maxit = Theta.maxit,
                                approx = FALSE, penalize.diagonal = TRUE, trace = FALSE)
  } else {
    Theta.out <- glasso::glasso(s = info$residual.cov, rho = lamTh, thr = Theta.thr, maxit = Theta.maxit,
                                approx = FALSE, penalize.diagonal = FALSE, trace = FALSE)
  }
  Theta <- (Theta.out$wi + t(Theta.out$wi))/2
  
  B.out <- updateBeta(Theta = Theta, B0 = info$B.init, lamB = lamB, eta = 0.8,
                      tolin = Beta.thr, maxitrin = Beta.maxit, info = info)
  
  if (verbose == 2) {
    cat("  lambda.Theta: ", lamTh, ";  lambda.Beta: ", lamB, "\n")
    cat("  Iterations for updating Theta ", Theta.out$niter, "\n")
    cat("  Iterations for updating Beta: ", B.out$it.final, "\n\n")
  }
  
  return(list(Beta = B.out$Bhat, Theta = Theta))
}


