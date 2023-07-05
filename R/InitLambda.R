InitLambda <- function(lamB, lamTh, n.tr, init.obj, 
                       n.lamB, n.lamTh,
                       lamB.min.ratio, lamTh.min.ratio,
                       lamB.scale.factor, lamTh.scale.factor) {
  if (!is.null(lamB) & !is.null(lamTh)) {
    lamB.vec <- sort(lamB, decreasing = TRUE)
    lamTh.vec <- sort(lamTh, decreasing = TRUE)
    lamB.vec.long <- rep(lamB.vec, each = length(lamTh.vec))
    lamTh.vec.long <- NULL
    for (l in 1:length(lamB.vec)) {
      if (l%%2 == 0) {
        lamTh.vec.long <- c(lamTh.vec.long, rev(lamTh.vec))
      } else {
        lamTh.vec.long <- c(lamTh.vec.long, lamTh.vec)
      }
    }
    
  } else {
    if (n.tr > nrow(init.obj$B.init)) {
      if (is.null(lamB.min.ratio)) { lamB.min.ratio <- 1e-4 }
      if (is.null(n.lamB)) { n.lamB <- 40 }
    } else {
      if (is.null(lamB.min.ratio)) { lamB.min.ratio <- 1e-3 }
      if (is.null(n.lamB)) { n.lamB <- 30 }
    }
    if (n.tr > ncol(init.obj$B.init)) {
      if (is.null(lamTh.min.ratio)) { lamTh.min.ratio <- 1e-4 }
      if (is.null(n.lamTh)) { n.lamTh <- 40 }
    } else {
      if (is.null(lamTh.min.ratio)) { lamTh.min.ratio <- 1e-3 }
      if (is.null(n.lamTh)) { n.lamTh <- 30 }
    }
    
    if (is.null(lamB) & is.null(lamTh)) {
      log.lamB.max <- log10(init.obj$lamB.max * lamB.scale.factor)
      log.lamB.min <- log10(init.obj$lamB.max * lamB.scale.factor * lamB.min.ratio)
      lamB.vec <- 10^(seq(from = log.lamB.max, to = log.lamB.min, length.out = n.lamB))
      lamB.vec.long <- rep(lamB.vec, each = n.lamTh)
      
      log.lamTh.max <- log10(init.obj$lamTh.max * lamTh.scale.factor)
      log.lamTh.min <- log10(init.obj$lamTh.max * lamTh.scale.factor * lamTh.min.ratio)
      lamTh.vec <- 10^(seq(from = log.lamTh.max, to = log.lamTh.min, length.out = n.lamTh))
      lamTh.vec.long <- NULL
      for (l in 1:n.lamB) {
        if (l%%2 == 0) {
          lamTh.vec.long <- c(lamTh.vec.long, rev(lamTh.vec))
        } else {
          lamTh.vec.long <- c(lamTh.vec.long, lamTh.vec)
        }
      }
      
    } else if (is.null(lamB) & !is.null(lamTh)) {
      log.lamB.max <- log10(init.obj$lamB.max * lamB.scale.factor)
      log.lamB.min <- log10(init.obj$lamB.max * lamB.scale.factor * lamB.min.ratio)
      lamB.vec <- 10^(seq(from = log.lamB.max, to = log.lamB.min, length.out = n.lamB))
      lamB.vec.long <- rep(lamB.vec, each = length(lamTh))
      
      lamTh.vec <- sort(lamTh, decreasing = TRUE)  ## lamTh.vec should be decreasing from largest to smallest
      lamTh.vec.long <- NULL
      for (l in 1:length(lamB.vec)) {
        if (l%%2 == 0) {
          lamTh.vec.long <- c(lamTh.vec.long, rev(lamTh.vec))
        } else {
          lamTh.vec.long <- c(lamTh.vec.long, lamTh.vec)
        }
      }
      
    } else if (!is.null(lamB) & is.null(lamTh)) {
      log.lamTh.max <- log10(init.obj$lamTh.max * lamTh.scale.factor)
      log.lamTh.min <- log10(init.obj$lamTh.max * lamTh.scale.factor * lamTh.min.ratio)
      lamTh.vec <- 10^(seq(from = log.lamTh.max, to = log.lamTh.min, length.out = n.lamTh))
      lamTh.vec.long <- NULL
      for (l in 1:length(lamB)) {
        if (l%%2 == 0) {
          lamTh.vec.long <- c(lamTh.vec.long, rev(lamTh.vec))
        } else {
          lamTh.vec.long <- c(lamTh.vec.long, lamTh.vec)
        }
      }
      
      lamB.vec <- sort(lamB, decreasing = TRUE)  ## lamB.vec should be decreasing from largest to smallest
      lamB.vec.long <- rep(lamB.vec, each = length(lamTh.vec))
    }
  }
  
  return(list(lamB.vec = lamB.vec.long, lamTh.vec = lamTh.vec.long))
}

