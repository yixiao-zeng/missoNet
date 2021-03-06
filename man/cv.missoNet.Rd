% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.missoNet.R
\name{cv.missoNet}
\alias{cv.missoNet}
\title{Cross-validation for missoNet}
\usage{
cv.missoNet(
  X,
  Y,
  kfold = 5,
  rho = NULL,
  lambda.Beta = NULL,
  lambda.Theta = NULL,
  lamBeta.min.ratio = 0.01,
  lamTheta.min.ratio = 0.01,
  lamBeta.scale.factor = 1,
  lamTheta.scale.factor = 1,
  n.lamBeta = round(-log10(lamBeta.min.ratio) * 15),
  n.lamTheta = round(-log10(lamTheta.min.ratio) * 10),
  Beta.maxit = 1000,
  Beta.thr = 1e-04,
  eta = 0.8,
  Theta.maxit = 1000,
  Theta.thr = 1e-04,
  eps = 1e-08,
  diag.penalty.factor = NULL,
  standardize = TRUE,
  standardize.response = TRUE,
  fit.1se = FALSE,
  fit.relax = FALSE,
  permute = FALSE,
  with.seed = NULL,
  parallel = FALSE,
  cpus,
  verbose = 1
)
}
\arguments{
\item{X}{Numeric predictor matrix (n by p): columns correspond to predictor variables and rows correspond to samples. Missing values are not allowed. Do not include a column of ones.}

\item{Y}{Numeric response matrix (n by q): columns correspond to response variables and rows correspond to samples. Missing values should be coded as \code{NA} or \code{NaN}.}

\item{kfold}{Number of folds -- default is 5.}

\item{rho}{(Optional) A scalar or a numeric vector of length q: a user supplied missing probability for response variables. Default is \code{rho = NULL} and the program will compute the empirical missing rates for columns of \code{Y} and use them as the missing probability.}

\item{lambda.Beta}{(Optional) Numeric vector: a user supplied regularization parameter sequence for \code{Beta} from which the CV procedure searches; default is \code{lambda.Beta = NULL} and the program computes its own lambda sequence based on \code{n.lamBeta} and \code{lamBeta.min.ratio}, supplying a vector overrides this. WARNING: avoid supplying a single value, missoNet relies on warm starts for acceleration and accuracy. The sequence supplied will be automatically arranged in descending order internally.}

\item{lambda.Theta}{(Optional) Numeric vector: a user supplied regularization parameter sequence for \code{Theta} from which the CV procedure searches; default is \code{lambda.Theta = NULL} and the program computes its own lambda sequence based on \code{n.lamTheta} and \code{lamTheta.min.ratio}, supplying a vector overrides this. WARNING: avoid supplying a single value, missoNet relies on warm starts for acceleration and accuracy. The sequence supplied will be automatically arranged in descending order internally.}

\item{lamBeta.min.ratio}{The smallest value for \code{lambda.Beta} will be data-derived \code{max(lambda.Beta)} multiplied by \code{lamBeta.min.ratio}. The default depends on the sample size n relative to the number of predictors p. If n > p, default is 0.01; if n <= p, default is 0.1. A very small value of \code{lamBeta.min.ratio} will lead to long computation time and may cause a saturated fit in the n < p case. Only needed when \code{lambda.Beta = NULL}.}

\item{lamTheta.min.ratio}{The smallest value for \code{lambda.Theta} will be data-derived \code{max(lambda.Theta)} multiplied by \code{lamTheta.min.ratio}; default is 0.01. A very small value of \code{lamTheta.min.ratio} will lead to long computation time. Only needed when \code{lambda.Theta = NULL}.}

\item{lamBeta.scale.factor}{A multiplication factor for scaling the entire \code{lambda.Beta} sequence; defaults is 1. A typical usage scenario is when the optimal lambda approaches the boundary of the search area. Only needed when \code{lambda.Beta = NULL}.}

\item{lamTheta.scale.factor}{A multiplication factor for scaling the entire \code{lambda.Theta} sequence; defaults is 1. A typical usage scenario is when the optimal lambda approaches the boundary of the search area. Only needed when \code{lambda.Theta = NULL}.}

\item{n.lamBeta}{The number of \code{lambda.Beta} values. By default, the program estimates a number based on \code{lamBeta.min.ratio}. An excessively large number results in a significant increase in computation time. Only needed when \code{lambda.Beta = NULL}.}

\item{n.lamTheta}{The number of \code{lambda.Theta} values. By default, the program estimates a number based on \code{lamTheta.min.ratio}. An excessively large number results in a significant increase in computation time. Only needed when \code{lambda.Theta = NULL}.}

\item{Beta.maxit}{The maximum number of internal iterations allowed for updating \code{Beta}. Default is \code{Beta.maxit = 1000}.}

\item{Beta.thr}{The convergence threshold for updating \code{Beta}; default is \code{Beta.thr = 1e-4}. Iterations stop when absolute parameter change is less than \code{Beta.thr * sum(abs(Beta))}.}

\item{eta}{Backtracking line search shrinkage factor, default is \code{eta = 0.8}. You may want to choose a more appropriate \code{eta} for a faster \code{Beta} convergence based on your dataset. Note that \code{eta} must be greater than 0 and smaller than 1.}

\item{Theta.maxit}{The maximum number of internal iterations allowed for updating \code{Theta}. Default is \code{Theta.maxit = 1000}.}

\item{Theta.thr}{The convergence threshold for updating \code{Theta}; default is \code{Theta.thr = 1e-4}. Iterations stop when average absolute parameter change is less than \code{Theta.thr * ave(abs(offdiag(Sigma)))}.}

\item{eps}{A numeric tolerance level for L1 projection; default is \code{eps = 1e-8}. If any of the eigenvalues is less than the given tolerance, the unbiased estimate of covariance is projected onto L1 ball to have \code{min(eigen(Sigma)$value) == eps}.}

\item{diag.penalty.factor}{Numeric: a separate penalty factor for the diagonal entries of \code{Theta}. This is a number that multiplies \code{lambda.Theta} to allow differential shrinkage. Default is \code{NULL} and the program computes it based on an initial estimate of \code{Theta}. Can be 0 for no shrinkage. Only needed when n <= p.}

\item{standardize}{Logical: should the columns of \code{X} be standardized so each has unit length and zero average; default is \code{TRUE}. The parameter estimates will be returned on the original scale. If \code{X} has been standardized prior to fitting the model, you might not wish to standardize.}

\item{standardize.response}{Logical: should the columns of \code{Y} be standardized so each has unit length and zero average; default is \code{TRUE}. The parameter estimates will be returned on the original scale. If \code{Y} has been standardized prior to fitting the model, you might not wish to standardize.}

\item{fit.1se}{Logical: default is \code{FALSE}. Should the model be refitted with the largest \code{lambda.Beta} according to the one-standard-error rule?}

\item{fit.relax}{Logical: default is \code{FALSE}. If \code{TRUE}, the program will re-estimate a relaxed graph (\code{Theta}) without penalization, which could be useful for network analysis. WARNING: there may be convergence issues if the residual covariance matrix is not of full rank.}

\item{permute}{Logical: should the subject indices for CV be permuted? Default is \code{FALSE}.}

\item{with.seed}{A random seed for permutation.}

\item{parallel}{Logical: default is \code{FALSE}. If \code{TRUE}, use parallel cluster to fit each fold.}

\item{cpus}{Number of cores for parallelization. Only needed when \code{parallel = TRUE}.}

\item{verbose}{Value of 0, 1 or 2. 0 -- silent; 1 -- limited tracing; 2 -- detailed tracing. Limited tracing if \code{parallel = TRUE}.}
}
\value{
This function returns a cv.missoNet object containing a \code{list} of estimates and working parameters:
\itemize{
\item \code{est.min}: estimates using the pair of regularization parameters that gives the minimum cross-validation error.
\item \code{est.1se}: estimates using the pair of regularization parameters such that error is within one-standard-error of the minimum, if \code{fit.1se = TRUE}.
\item \code{lambda.Beta.vec}: a long vector of length \code{n.lamBeta * n.lamTheta} from which the CV procedure searches.
\item \code{lambda.Theta.vec}: a long vector of length \code{n.lamBeta * n.lamTheta} from which the CV procedure searches.
\item \code{relax.graph}: an estimate of network without penalization if \code{fit.relax = TRUE}.
\item \code{rho}: a vector of length q: working missing probability.
\item \code{fold.index}: subject indices for cross-validation split.
\item \code{cvm}: mean cross-validated error that has a one to one correspondence with \code{lambda.Beta.vec} and \code{lambda.Theta.vec}.
\item \code{cvup}: upper cross-validated error.
\item \code{cvlo}: lower cross-validated error.
\item \code{penalize.diagonal}: logical: whether the diagonal of \code{Theta} was penalized.
\item \code{diag.penalty.factor}: penalty factor for the diagonal when \code{penalize.diagonal} returned as \code{TRUE}.
}
}
\description{
Cross-validation for missoNet
}
