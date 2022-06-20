plot.heatmap <- function(cv.missoNet.obj) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    if (!require("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("ComplexHeatmap")
    requireNamespace("ComplexHeatmap", quietly = TRUE)
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    install.packages("circlize")
    requireNamespace("circlize", quietly = TRUE)
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    install.packages("grid")
    requireNamespace("grid", quietly = TRUE)
  }
  
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
  rownames(cvm) <- sprintf("%.3f", lamTh.vec)
  colnames(cvm) <- sprintf("%.3f", lamB.vec)
  
  ## Flip vertically
  cvm <- apply(cvm, 1, rev)
  ## Transpose
  cvm <- t(cvm)
  
  col <- circlize::colorRamp2(quantile(cv.missoNet.obj$cvm, c(0, 0.2, 0.4, 0.6, 0.8, 1)), c("blue", "cyan", "green", "yellow", "orange", "red"))
  ComplexHeatmap::Heatmap(cvm, col = col, border = TRUE,
                          row_title = expression(lambda[Theta]), row_title_side = "right", row_title_rot = 0,
                          row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9),
                          column_title = expression(lambda[Beta]), column_names_gp = grid::gpar(fontsize = 9),
                          cluster_rows = FALSE, cluster_columns = FALSE,
                          name = "CV.Error", # title of legend
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            if(j == which(sort(lamB.vec, decreasing=FALSE) == cv.missoNet.obj$est.min$lambda.Beta) 
                               & i == which(sort(lamTh.vec, decreasing=TRUE) == cv.missoNet.obj$est.min$lambda.Theta)) {
                              grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = "white", fill = NA))
                            }
                          }
  )
  
}


