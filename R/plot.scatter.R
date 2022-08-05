plot.scatter <- function(cv.missoNet.obj, detailed.axes) {
  if (!requireNamespace("plot3D", quietly = TRUE)) {
    install.packages("plot3D")
    requireNamespace("plot3D", quietly = TRUE)
  }
  
  log.lamB.vec <- log10(cv.missoNet.obj$lambda.Beta.vec)
  log.lamTh.vec <- log10(cv.missoNet.obj$lambda.Theta.vec)
  
  old <- par(cex = 0.9, mai = c(0.6, 0.2, 0.3, 0.5))
  on.exit(par(old), add = TRUE)
  if (detailed.axes) {
    plot3D::scatter3D(x = log.lamB.vec, y = log.lamTh.vec, z = cv.missoNet.obj$cvm,
                      pch = 16, cex = 0.83, col = plot3D::ramp.col(col=c("blue","cyan","green","yellow","orange","red"), alpha=0.8),
                      xlab = "log10(lambda.Beta)", ylab = "log10(lambda.Theta)", zlab = "CV.Error", clab = c("Magnitude"),
                      ticktype = "detailed", theta = -40, d = 5, phi = 20,
                      colkey = list(length = 0.63, width = 0.66, cex.clab = 0.93)
    )} else {
      plot3D::scatter3D(x = log.lamB.vec, y = log.lamTh.vec, z = cv.missoNet.obj$cvm,
                        pch = 16, cex = 0.83, col = plot3D::ramp.col(col=c("blue","cyan","green","yellow","orange","red"), alpha=0.8),
                        xlab = "log10(lambda.Beta)", ylab = "log10(lambda.Theta)", zlab = "CV.Error", clab = c("Magnitude"),
                        theta = -40, d = 5, phi = 20,
                        colkey = list(length = 0.63, width = 0.66, cex.clab = 0.93)
      )}
}

