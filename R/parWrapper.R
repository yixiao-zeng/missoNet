parWrapper <- function(k, X, Y, init.obj, rho, ind, kfold, lamTh.vec, lamB.vec,
                       Beta.maxit, Beta.thr, Theta.maxit, Theta.thr, eps, eta) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  err.fold <- rep(0, length(lamTh.vec))
  Beta.warm.fold <- lapply(1:length(lamTh.vec), function(x){matrix(0,p,q)})
  
  foldind <- ind[(1+floor((k-1)*n/kfold)) : floor(k*n/kfold)]
  X.tr <- X[-foldind, ]
  Y.tr <- Y[-foldind, ]
  X.va <- X[foldind, , drop = FALSE]
  Y.va <- Y[foldind, , drop = FALSE]
  n.tr <- nrow(Y.tr)
  
  if (is.null(rho)) {
    rho.vec <- apply(Y.tr, 2, function(x) {
      sum(is.na(x))/n.tr
    })
  } else { rho.vec <- init.obj$rho.vec }
  rho.mat.1 <- t(matrix(rep(1 - rho.vec, p), q, p))  # pxq
  rho.mat.2 <- matrix(1 - rho.vec, q, 1) %*% matrix(1 - rho.vec, 1, q)
  diag(rho.mat.2) <- 1 - rho.vec  # qxq
  
  mx.tr <- apply(X.tr, 2, mean)
  X.tr <- scale(X.tr, center = mx.tr, scale = init.obj$sdx)
  X.va <- scale(X.va, center = mx.tr, scale = init.obj$sdx)
  
  my.tr <- apply(Y.tr, 2, mean, na.rm = TRUE)
  Y.tr <- scale(Y.tr, center = my.tr, scale = init.obj$sdy)
  Y.va <- scale(Y.va, center = my.tr, scale = init.obj$sdy)
  Z.tr <- Y.tr
  Z.tr[is.na(Z.tr)] <- 0
  
  info <- NULL
  info$n <- n.tr
  info$q <- q
  info$penalize.diagonal <- init.obj$penalize.diagonal
  info$xtx <- crossprod(X.tr)
  info$til.xty <- crossprod(X.tr, Z.tr)/rho.mat.1
  
  info.update <- NULL
  info.update$B.init <- init.obj$B.init * init.obj$sdx
  Beta.thr.rescale <- Beta.thr * sum(abs(info.update$B.init))
  E.tr <- Y.tr - X.tr %*% info.update$B.init
  info.update$residual.cov <- getResidual(E = E.tr, n = n.tr, rho.mat = rho.mat.2, eps = eps)
  
  for (i in 1:length(lamTh.vec)) {
    info.update$B.init <- update.missoNet(lamTh = lamTh.vec[i], lamB = lamB.vec[i],
                                          Beta.maxit = Beta.maxit, Beta.thr = Beta.thr.rescale,
                                          Theta.maxit = Theta.maxit, Theta.thr = Theta.thr,
                                          verbose = 0, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                                          info = info, info.update = info.update, under.cv = TRUE)
    Beta.thr.rescale <- Beta.thr * sum(abs(info.update$B.init))
    E.tr <- Y.tr - X.tr %*% info.update$B.init
    info.update$residual.cov <- getResidual(E = E.tr, n = n.tr, rho.mat = rho.mat.2, eps = eps)
    
    E.va.sq <- (Y.va - X.va %*% info.update$B.init)^2
    err.fold[i] <- mean(E.va.sq, na.rm = TRUE)
    Beta.warm.fold[[i]] <- info.update$B.init
  }
  
  return(list(err.fold = err.fold, Beta.warm.fold = Beta.warm.fold))
}

