plot.scatter <- function(cv.missoNet.obj, detailed.axes, plt.surf, ...) {
  col <- circlize::colorRamp2(quantile(cv.missoNet.obj$cvm, c(0, 0.2, 0.4, 0.6, 0.8, 1)), c("blue", "cyan", "green", "yellow", "orange", "red"))
  if (isTRUE(detailed.axes)) {
    lab <- c(10, 10, 10)
  } else { lab <- par("lab") }
  
  s3d <- scatterplot3d::scatterplot3d(log10(cv.missoNet.obj$lambda.Beta.vec), log10(cv.missoNet.obj$lambda.Theta.vec), cv.missoNet.obj$cvm,
                                      type = "p", pch = 20, cex.symbols = 0.9, color = col(cv.missoNet.obj$cvm), lab = lab,
                                      zlab = "CV.Error", xlab = expression(log10(lambda[Beta])), ylab = expression(log10(lambda[Theta])), ...)
  
  if (isTRUE(plt.surf)) {
    lamB.vec <- sort(unique(log10(cv.missoNet.obj$lambda.Beta.vec)), decreasing = FALSE)
    lamTh.vec <- sort(unique(log10(cv.missoNet.obj$lambda.Theta.vec)), decreasing = FALSE)
    cvm <- NULL
    for (l in 1:length(lamB.vec)) {
      if (l%%2 == 0) {
        cvm <- cbind(cvm, rev(cv.missoNet.obj$cvm[((l - 1) * length(lamTh.vec) + 1):(l * length(lamTh.vec))]))
      } else {
        cvm <- cbind(cvm, cv.missoNet.obj$cvm[((l - 1) * length(lamTh.vec) + 1):(l * length(lamTh.vec))])
      }
    }
    cvm <- apply(cvm, 1, rev)
    cvm <- t(cvm)
    
    for (i in 1:length(lamB.vec))
      s3d$points3d(rep(lamB.vec[i], length(lamTh.vec)), rev(lamTh.vec), cvm[ ,i], type = "l", lty = "dotted")
    for(i in length(lamTh.vec):1)
      s3d$points3d(lamB.vec, rep(lamTh.vec[i], length(lamB.vec)), cvm[length(lamTh.vec)-i+1, ], type = "l", lty = "dotted")
  }
}
