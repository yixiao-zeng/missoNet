plot.heatmap <- function(cv.missoNet.obj, detailed.axes, ...) {
  lamB.vec <- sort(unique(cv.missoNet.obj$lambda.Beta.vec), decreasing = TRUE)
  lamTh.vec <- sort(unique(cv.missoNet.obj$lambda.Theta.vec), decreasing = TRUE)
  cvm <- NULL
  for (l in 1:length(lamB.vec)) {
    if (l%%2 == 0) {
      cvm <- cbind(cvm, rev(cv.missoNet.obj$cvm[((l - 1) * length(lamTh.vec) + 1):(l * length(lamTh.vec))]))
    } else {
      cvm <- cbind(cvm, cv.missoNet.obj$cvm[((l - 1) * length(lamTh.vec) + 1):(l * length(lamTh.vec))])
    }
  }
  if (detailed.axes) {
    rownames(cvm) <- sprintf("%.3f", lamTh.vec)
    colnames(cvm) <- sprintf("%.3f", lamB.vec)
  } else {
    rownames(cvm) <- rep(" ", length(lamTh.vec))
    colnames(cvm) <- rep(" ", length(lamB.vec))
    pos <- c(1, floor(length(lamTh.vec) * 0.1), floor(length(lamTh.vec) * 0.2),
             floor(length(lamTh.vec) * 0.3), floor(length(lamTh.vec) * 0.4),
             floor(length(lamTh.vec) * 0.5), floor(length(lamTh.vec) * 0.6),
             floor(length(lamTh.vec) * 0.7), floor(length(lamTh.vec) * 0.8),
             floor(length(lamTh.vec) * 0.9), length(lamTh.vec))
    rownames(cvm)[pos] <- sprintf("%.3f", lamTh.vec[pos])
    pos <- c(1, floor(length(lamB.vec) * 0.1), floor(length(lamB.vec) * 0.2),
             floor(length(lamB.vec) * 0.3), floor(length(lamB.vec) * 0.4),
             floor(length(lamB.vec) * 0.5), floor(length(lamB.vec) * 0.6),
             floor(length(lamB.vec) * 0.7), floor(length(lamB.vec) * 0.8),
             floor(length(lamB.vec) * 0.9), length(lamB.vec))
    colnames(cvm)[pos] <- sprintf("%.3f", lamB.vec[pos])
  }
  
  ## Flip vertically
  cvm <- apply(cvm, 1, rev)
  ## Transpose
  cvm <- t(cvm)
  
  col <- circlize::colorRamp2(quantile(cv.missoNet.obj$cvm, c(0, 0.2, 0.4, 0.6, 0.8, 1)), c("blue", "cyan", "green", "yellow", "orange", "red"))
  ComplexHeatmap::Heatmap(cvm, col = col, border = TRUE,
                          row_title = expression(lambda[Theta]), row_title_side = "right", row_title_rot = 0,
                          row_names_side = "left", row_names_gp = gpar(fontsize = 8.4),
                          column_title = expression(lambda[Beta]), column_names_gp = gpar(fontsize = 8.4),
                          cluster_rows = FALSE, cluster_columns = FALSE,
                          name = "CV.Error", # title of legend
                          ...,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            if (j == which(sort(lamB.vec, decreasing = FALSE) == cv.missoNet.obj$est.min$lambda.Beta) 
                                & i == which(sort(lamTh.vec, decreasing = TRUE) == cv.missoNet.obj$est.min$lambda.Theta)) {
                              grid.rect(x = x, y = y, width = width, height = height, gp = gpar(lwd = 1.5, col = "white", fill = NA))
                            }
                            if (!is.null(cv.missoNet.obj$est.1se.B)) {
                              if (j == which(sort(lamB.vec, decreasing = FALSE) == cv.missoNet.obj$est.1se.B$lambda.Beta) 
                                  & i == which(sort(lamTh.vec, decreasing = TRUE) == cv.missoNet.obj$est.1se.B$lambda.Theta)) {
                                grid.rect(x = x, y = y, width = width, height = height, gp = gpar(lwd = 1.5, lty = "dashed", col = "white", fill = NA))
                              }
                            }
                            if (!is.null(cv.missoNet.obj$est.1se.Tht)) {
                              if (j == which(sort(lamB.vec, decreasing = FALSE) == cv.missoNet.obj$est.1se.Tht$lambda.Beta) 
                                  & i == which(sort(lamTh.vec, decreasing = TRUE) == cv.missoNet.obj$est.1se.Tht$lambda.Theta)) {
                                grid.rect(x = x, y = y, width = width, height = height, gp = gpar(lwd = 1.5, lty = "dashed", col = "white", fill = NA))
                              }
                            }
                          }
  )
}
