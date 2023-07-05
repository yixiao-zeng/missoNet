boundaryCheck <- function(lambda.Theta, lambda.Beta, lamTh.vec, lamB.vec, lamTh.min, lamB.min, margin) {
  bound <- quantile(unique(lamTh.vec), c(margin, 1 - margin))
  if (lamTh.min >= bound[2]) {
    if (is.null(lambda.Theta)) {
      warning("\nThe optimal `lambda.Theta` is close to the upper boundary of the search scope, 
try to provide a larger value for the following argument:\n
    1. lamTheta.scale.factor\n")
    } else {
      warning("\nThe optimal `lambda.Theta` is close to the upper boundary of the search scope, 
try to provide a new sequence covering larger values for the following argument:\n
    1. lambda.Theta\n")
    }
  } else if (lamTh.min <= bound[1]) {
    if (is.null(lambda.Theta)) {
      warning("\nThe optimal `lambda.Theta` is close to the lower boundary of the search scope, 
try to provide a smaller value for one of the following arguments (or both):\n
    1. lamTheta.scale.factor (> 0)
    2. lamTheta.min.ratio (> 0)\n")
    } else {
      warning("\nThe optimal `lambda.Theta` is close to the lower boundary of the search scope, 
try to provide a new sequence covering smaller values for the following argument:\n
    1. lambda.Theta\n")
    }
  }
  
  bound <- quantile(unique(lamB.vec), c(margin, 1 - margin))
  if (lamB.min >= bound[2]) {
    if (is.null(lambda.Beta)) {
      warning("\nThe optimal `lambda.Beta` is close to the upper boundary of the search scope, 
try to provide a larger value for the following argument:\n
      1. lamBeta.scale.factor\n")
    } else {
      warning("\nThe optimal `lambda.Beta` is close to the upper boundary of the search scope, 
try to provide a new sequence covering larger values for the following argument:\n
      1. lambda.Beta\n")
    }
  } else if (lamB.min <= bound[1]) {
    if (is.null(lambda.Beta)) {
      warning("\nThe optimal `lambda.Beta` is close to the lower boundary of the search scope, 
try to provide a smaller value for one of the following arguments (or both):\n
      1. lamBeta.scale.factor (> 0)
      2. lamBeta.min.ratio (> 0)\n")
    } else {
      warning("\nThe optimal `lambda.Beta` is close to the lower boundary of the search scope, 
try to provide a new sequence covering smaller values for the following argument:\n
      1. lambda.Beta\n")
    }
  }
}

