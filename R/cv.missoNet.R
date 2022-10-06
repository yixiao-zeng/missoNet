#' Cross-validation for missoNet
#'
#' This function performs k-fold cross-validation for \sQuote{\code{\link{missoNet}}}. The regularization path is computed 
#' for all possible combinations of values given in the two regularization parameter sequences, namely \eqn{\lambda_B} and \eqn{\lambda_\Theta}. 
#' \sQuote{\code{cv.missoNet}} will select the most suitable model among all cross-validated fits along the path.
#' See the details of \sQuote{\code{\link{missoNet}}} for the model definition. 
#' To help users, the \sQuote{\code{cv.missoNet}} function is designed to automatically determine the likely ranges 
#' of the regularization parameters over which the cross-validation searches.
#' 
#' The \sQuote{\code{cv.missoNet}} function fits \sQuote{\code{\link{missoNet}}} models (\code{'kfold'} \code{*} \code{'n.lamBeta'} \code{*} \code{'n.lamTheta'}) 
#' times in the whole cross-validation process. That is, for the \eqn{k}th-fold (\eqn{k=1,...,K}) computation, the models are fitted at each of 
#' the all (\code{'n.lamBeta'} \code{*} \code{'n.lamTheta'}) possible combinations of the regularization parameters (\eqn{\lambda_B}, \eqn{\lambda_\Theta}), with the \eqn{k}th 
#' fold of the training data omitted. The errors are accumulated, and the averaged errors as well as the standard deviations are computed over all folds. 
#' Note that the results of \sQuote{\code{cv.missoNet}} are random, since the samples are randomly split into k-folds. Users can eliminate this randomness 
#' by setting \code{'permute = FALSE'}, or by explicitly assigning a seed to the permutation of samples.
#' 
#' A user-supplied sequence for \{\eqn{\lambda_B}\} and/or \{\eqn{\lambda_\Theta}\} is permitted, 
#' otherwise the program computes an appropriate range of values for the regularization parameters using other control arguments.
#' Note that \sQuote{\code{cv.missoNet}} standardizes \code{'X'} and \code{'Y'} to have unit variances before computing its \eqn{\lambda} 
#' sequences (and then unstandardizes the resulting coefficients); if you wish to reproduce/compare results with those of other softwares, 
#' it is best to supply pre-standardized \code{'X'} and \code{'Y'}. If the algorithm is running slowly, track its progress with \code{'verbose = 2'}. 
#' In most cases, choosing a sparser grid for the regularization parameters (e.g. smaller \code{'n.lamBeta'} and/or \code{'n.lamTheta'}) or setting \code{'Beta.thr = 1.0E-3'} (or even bigger) 
#' allows the algorithm to make faster progress.
#' 
#' After cross-validation, the regression coefficient matrix \eqn{\mathbf{B}} and the precision matrix \eqn{\mathbf{\Theta}} can be 
#' estimated at three special \eqn{\lambda} pairs, by reapplying \sQuote{\code{missoNet}} to the entire training dataset:
#' \enumerate{
#'   \item "\code{lambda.min}" := (\eqn{{\lambda_B}_\mathrm{min}, {\lambda_\Theta}_\mathrm{min}}), at which the minimum mean cross-validated error is achieved;
#'   \item "\code{lambda.1se.Beta}" := (\eqn{{\lambda_B}_\mathrm{1se}, {\lambda_\Theta}_\mathrm{min}}), where \eqn{{\lambda_B}_\mathrm{1se}} is the largest \eqn{\lambda_B} at which the mean cross-validated error is within one standard error of the minimum;
#'   \item "\code{lambda.1se.Theta}" := (\eqn{{\lambda_B}_\mathrm{min}, {\lambda_\Theta}_\mathrm{1se}}), where \eqn{{\lambda_\Theta}_\mathrm{1se}} is the largest \eqn{\lambda_\Theta} at which the mean cross-validated error is within one standard error of the minimum.
#' }
#' The corresponding estimates, along with the \eqn{\lambda} values, are stored in three separate lists inside the returned object: 
#' \code{'est.min'}, \code{'est.1se.B'} and \code{'est.1se.Tht'} (\code{'est.1se.B'} and \code{'est.1se.Tht'} are \code{'NULL'} 
#' unless the argument \code{'fit.1se = TRUE'} when calling \sQuote{\code{cv.missoNet}}).
#' 
#' The \sQuote{\code{cv.missoNet}} function returns an R object of S3 class \code{'cv.missoNet'} for which there are a set of accessory functions available.
#' The plotting function \sQuote{\code{\link{plot.cv.missoNet}}} can be used to graphically identify the optimal pair of the regularization parameters, 
#' and the prediction function \sQuote{\code{\link{predict.cv.missoNet}}} can be used to make predictions of response values given new input \code{'X'}. 
#' See the vignette for examples.
#' 
#' @param X Numeric predictor matrix (\eqn{n\times p}): columns correspond to predictor variables and rows correspond to samples. Missing values are not allowed. There is no need for centering or scaling of the variables. \code{'X'} should not include a column of ones for an intercept.
#' @param Y Numeric response matrix (\eqn{n\times q}): columns correspond to response variables and rows correspond to samples. Missing values should be coded as either \code{'NA'}s or \code{'NaN'}s. There is no need for centering or scaling of the variables.
#' @param kfold Number of folds for cross-validation -- the default is \code{'5'}.
#' @param rho (Optional) A scalar or a numeric vector of length \eqn{q}: the elements are user-supplied probabilities of missingness for the response variables. The default is \code{'rho = NULL'} and the program will compute the empirical missing rates for each of the columns of \code{'Y'} and use them as the working missing probabilities. The default setting should suffice in most cases; misspecified missing probabilities would introduce biases into the model.
#' @param lambda.Beta (Optional) Numeric vector: a user-supplied sequence of non-negative values for \{\eqn{\lambda_B}\} penalizing the elements of the coefficient matrix \eqn{\mathbf{B}} among which the cross-validation procedure searches. The default is \code{'lambda.Beta = NULL'}, in which case the program computes an appropriate range of \eqn{\lambda_B} values using \code{'n.lamBeta'} and \code{'lamBeta.min.ratio'}. Supplying a vector overrides this default. Note that the supplied sequence will be automatically arranged, internally, in a descending order.
#' @param lambda.Theta (Optional) Numeric vector: a user-supplied sequence of non-negative values for \{\eqn{\lambda_\Theta}\} penalizing the (off-diagonal) elements of the precision matrix \eqn{\mathbf{\Theta}} among which the cross-validation procedure searches. The default is \code{'lambda.Theta = NULL'}, in which case the program computes an appropriate range of \eqn{\lambda_\Theta} values using \code{'n.lamTheta'} and \code{'lamTheta.min.ratio'}. Supplying a vector overrides this default. Note that the supplied sequence will be automatically arranged, internally, in a descending order. 
#' @param lamBeta.min.ratio The smallest value of \eqn{\lambda_B} is calculated as the data-derived \eqn{\mathrm{max}(\lambda_B)} multiplied by \code{'lamBeta.min.ratio'}. The default depends on the sample size, \eqn{n}, relative to the number of predictors, \eqn{p}. If \eqn{n > p}, the default is \code{'1.0E-4'}, otherwise it is \code{'1.0E-2'}. A very small value of \code{'lamBeta.min.ratio'} may significantly increase runtime and lead to a saturated fit in the \eqn{n \leq p} case. This is only needed when \code{'lambda.Beta = NULL'}.
#' @param lamTheta.min.ratio The smallest value of \eqn{\lambda_\Theta} is calculated as the data-derived \eqn{\mathrm{max}(\lambda_\Theta)} multiplied by \code{'lamTheta.min.ratio'}. The default depends on the sample size, \eqn{n}, relative to the number of responses, \eqn{q}. If \eqn{n > q}, the default is \code{'1.0E-4'}, otherwise it is \code{'1.0E-2'}. A very small value of \code{'lamTheta.min.ratio'} may significantly increase runtime and lead to a saturated fit in the \eqn{n \leq q} case. This is only needed when \code{'lambda.Theta = NULL'}.
#' @param n.lamBeta The number of \eqn{\lambda_B} values. If \eqn{n > p}, the default is \code{'40'}, otherwise it is \code{'20'}. Avoid supplying an excessively large number since the program will fit (\code{'n.lamBeta'} \code{*} \code{'n.lamTheta'}) models in total for each fold of the cross-validation. Typically we suggest \code{'n.lamBeta' = -log10('lamBeta.min.ratio') * c}, where \code{c} \eqn{\in} [\code{10}, \code{20}]. This is only needed when \code{'lambda.Beta = NULL'}.
#' @param n.lamTheta The number of \eqn{\lambda_\Theta} values. If \eqn{n > q}, the default is \code{'40'}, otherwise it is \code{'20'}. Avoid supplying an excessively large number since the program will fit (\code{'n.lamBeta'} \code{*} \code{'n.lamTheta'}) models in total for each fold of the cross-validation. Typically we suggest \code{'n.lamTheta' = -log10('lamTheta.min.ratio') * c}, where \code{c} \eqn{\in} [\code{10}, \code{20}]. This is only needed when \code{'lambda.Theta = NULL'}.
#' @param lamBeta.scale.factor A positive multiplication factor for scaling the entire \eqn{\lambda_B} sequence; the default is \code{'1'}. A typical usage is when the magnitudes of the auto-computed \eqn{\lambda_B} values are inappropriate. For example, this factor would be needed if the optimal value of \eqn{\lambda_B} selected by the cross-validation (i.e. \eqn{{\lambda_B}_\mathrm{min}} with the minimum cross-validated error) approaches either boundary of the search range. This is only needed when \code{'lambda.Beta = NULL'}.
#' @param lamTheta.scale.factor A positive multiplication factor for scaling the entire \eqn{\lambda_\Theta} sequence; the default is \code{'1'}. A typical usage is when the magnitudes of the auto-computed \eqn{\lambda_\Theta} values are inappropriate. For example, this factor would be needed if the optimal value of \eqn{\lambda_\Theta} selected by the cross-validation (i.e. \eqn{{\lambda_\Theta}_\mathrm{min}} with the minimum cross-validated error) approaches either boundary of the search range. This is only needed when \code{'lambda.Theta = NULL'}.
#' @param Beta.maxit The maximum number of iterations of the fast iterative shrinkage-thresholding algorithm (FISTA) for updating \eqn{\hat{\mathbf{B}}}. The default is \code{'Beta.maxit = 1000'}.
#' @param Beta.thr The convergence threshold of the FISTA algorithm for updating \eqn{\hat{\mathbf{B}}}; the default is \code{'Beta.thr = 1.0E-4'}. Iterations stop when the absolute parameter change is less than (\code{'Beta.thr'} \code{*} \code{sum(abs(}\eqn{\hat{\mathbf{B}}}\code{))}).
#' @param eta The backtracking line search shrinkage factor; the default is \code{'eta = 0.8'}. Most users will be able to use the default value, some experienced users may want to tune \code{'eta'} according to the properties of a specific dataset for a faster convergence of the FISTA algorithm. Note that \code{'eta'} must be in (0, 1).
#' @param Theta.maxit The maximum number of iterations of the \sQuote{\code{\link{glasso}}} algorithm for updating \eqn{\hat{\mathbf{\Theta}}}. The default is \code{'Theta.maxit = 1000'}.
#' @param Theta.thr The convergence threshold of the \sQuote{\code{\link{glasso}}} algorithm for updating \eqn{\hat{\mathbf{\Theta}}}; the default is \code{'Theta.thr = 1.0E-4'}. Iterations stop when the average absolute parameter change is less than (\code{'Theta.thr'} \code{*} \code{ave(abs(offdiag(}\eqn{\hat{\mathbf{\Sigma}}}\code{)))}), where \eqn{\hat{\mathbf{\Sigma}}} denotes the empirical working covariance matrix.
#' @param eps A numeric tolerance level for the L1 projection of the empirical covariance matrix; the default is \code{'eps = 1.0E-8'}. The empirical covariance matrix will be projected onto a L1 ball to have \code{min(eigen(}\eqn{\hat{\mathbf{\Sigma}}}\code{)$value)} == \code{'eps'}, if any of the eigenvalues is less than the specified tolerance. Most users will be able to use the default value.
#' @param penalize.diagonal Logical: should the diagonal elements of \eqn{\mathbf{\Theta}} be penalized? The default depends on the sample size, \eqn{n}, relative to the number of predictors and responses. If \eqn{n > \mathrm{max}(p, q)}, the default is \code{'TRUE'}, otherwise it is set to \code{'FALSE'}. Most users will be able to use the default setting.
#' @param diag.penalty.factor Numeric: a separate penalty multiplication factor for the diagonal elements of \eqn{\mathbf{\Theta}} when \code{'penalize.diagonal = TRUE'}. \eqn{\lambda_\Theta} is multiplied by this number to allow a differential shrinkage of the diagonal elements. The default is \code{'NULL'} and the program will guess a value based on an initial estimate of \eqn{\mathbf{\Theta}}. This factor could be \code{'0'} for no shrinkage (equivalent to \code{'penalize.diagonal = FALSE'}). Most users will be able to use the default value.
#' @param standardize Logical: should the columns of \code{'X'} be standardized so each has unit variance? The default is \code{'TRUE'}. The estimated results will always be returned on the original scale. \sQuote{\code{cv.missoNet}} computes appropriate \eqn{\lambda} sequences relying on standardization, if \code{'X'} has been standardized prior to fitting the model, you might not wish to standardize it inside the algorithm.
#' @param standardize.response Logical: should the columns of \code{'Y'} be standardized so each has unit variance? The default is \code{'TRUE'}. The estimated results will always be returned on the original scale. \sQuote{\code{cv.missoNet}} computes appropriate \eqn{\lambda} sequences relying on standardization, if \code{'Y'} has been standardized prior to fitting the model, you might not wish to standardize it inside the algorithm.
#' @param fit.1se Logical: the default is \code{'FALSE'}. If \code{'TRUE'}, two additional models will be fitted with the largest values of \eqn{\lambda_B} and \eqn{\lambda_\Theta} respectively at which the cross-validated error is within one standard error of the minimum.
#' @param fit.relax Logical: the default is \code{'FALSE'}. If \code{'TRUE'}, the program will re-estimate the edges in the active set (i.e. nonzero off-diagonal elements) of the network structure \eqn{\hat{\mathbf{\Theta}}} without penalization (\eqn{\lambda_\Theta=0}). This debiased estimate of \eqn{\mathbf{\Theta}} could be useful for some interdependency analyses. WARNING: there may be convergence issues if the empirical covariance matrix is not of full rank (e.g. \eqn{n < q)}).
#' @param permute Logical: should the subject indices for the cross-validation be permuted? The default is \code{'TRUE'}.
#' @param with.seed A random number seed for the permutation.
#' @param parallel Logical: the default is \code{'FALSE'}. If \code{'TRUE'}, the program uses clusters to compute the cross-validation folds in parallel. Must register parallel clusters beforehand, see examples below.
#' @param cl A cluster object created by \sQuote{\code{parallel::makeCluster}} for parallel evaluations. This is only needed when \code{'parallel = TRUE'}.
#' @param verbose Value of \code{'0'}, \code{'1'} or \code{'2'}. \code{'verbose = 0'} -- silent; \code{'verbose = 1'} (the default) -- limited tracing with progress bars; \code{'verbose = 2'} -- detailed tracing. Note that displaying the progress bars slightly increases the computation overhead compared to the silent mode. The detailed tracing should be used for monitoring progress only when the program runs extremely slowly, and it is not supported under \code{'parallel = TRUE'}.
#'
#' @return This function returns a \code{'cv.missoNet'} object containing a named \code{'list'} with all the ingredients of the cross-validated fit:
#' \item{\code{est.min}}{A \code{'list'} of results estimated at "\code{lambda.min}" := (\eqn{{\lambda_B}_\mathrm{min}, {\lambda_\Theta}_\mathrm{min}}) that gives the minimum mean cross-validated error. It consists of the following components:
#'   \itemize{
#'       \item \code{Beta}: the penalized estimate of the regression coefficient matrix \eqn{\hat{\mathbf{B}}} (\eqn{p\times q}).
#'       \item \code{Theta}: the penalized estimate of the precision matrix \eqn{\hat{\mathbf{\Theta}}} (\eqn{q\times q}).
#'       \item \code{mu}: a vector of length \eqn{q} storing the model intercept \eqn{\hat{\mu}}.
#'       \item \code{lambda.Beta}: the value of \eqn{\lambda_B} (i.e. \eqn{{\lambda_B}_\mathrm{min}}) used to fit the model.
#'       \item \code{lambda.Theta}: the value of \eqn{\lambda_\Theta} (i.e. \eqn{{\lambda_\Theta}_\mathrm{min}}) used to fit the model.
#'       \item \code{relax.net}: the relaxed (debiased) estimate of the conditional network structure \eqn{\hat{\mathbf{\Theta}}_\mathrm{rlx}} (\eqn{q\times q}) if \code{'fit.relax = TRUE'} when calling \sQuote{\code{cv.missoNet}}.
#'   }
#' }
#' \item{\code{est.1se.B}}{A \code{'list'} of results estimated at "\code{lambda.1se.Beta}" := (\eqn{{\lambda_B}_\mathrm{1se}, {\lambda_\Theta}_\mathrm{min}}) if \code{'fit.1se = TRUE'} when calling \sQuote{\code{cv.missoNet}}. "\code{lambda.1se.Beta}" refers to the largest \eqn{\lambda_B} at which the mean cross-validated error is within one standard error of the minimum, by fixing \eqn{\lambda_\Theta} at \eqn{{\lambda_\Theta}_\mathrm{min}}. This \code{'list'} consists of the same components as \code{'est.min'}.}
#' \item{\code{est.1se.Tht}}{A \code{'list'} of results estimated at "\code{lambda.1se.Theta}" := (\eqn{{\lambda_B}_\mathrm{min}, {\lambda_\Theta}_\mathrm{1se}}) if \code{'fit.1se = TRUE'} when calling \sQuote{\code{cv.missoNet}}. "\code{lambda.1se.Theta}" refers to the largest \eqn{\lambda_\Theta} at which the mean cross-validated error is within one standard error of the minimum, by fixing \eqn{\lambda_B} at \eqn{{\lambda_B}_\mathrm{min}}. This \code{'list'} consists of the same components as \code{'est.min'}.}
#' \item{\code{rho}}{A vector of length \eqn{q} storing the working missing probabilities for the \eqn{q} response variables.}
#' \item{\code{fold.index}}{The subject indices identifying which fold each observation is in.}
#' \item{\code{lambda.Beta.vec}}{A flattened vector of length (\code{'n.lamBeta'} \code{*} \code{'n.lamTheta'}) storing the \eqn{\lambda_B} values along the regularization path. More specifically, \code{'lambda.Beta.vec' = rep('lambda.Beta', each = 'n.lamTheta')}.}
#' \item{\code{lambda.Theta.vec}}{A flattened vector of length (\code{'n.lamBeta'} \code{*} \code{'n.lamTheta'}) storing the \eqn{\lambda_\Theta} values along the regularization path. More specifically, \code{'lambda.Theta.vec' = rep('lambda.Theta', times = 'n.lamBeta')}.}
#' \item{\code{cvm}}{A flattened vector of length (\code{'n.lamBeta'} \code{*} \code{'n.lamTheta'}) storing the (standardized) mean cross-validated errors along the regularization path.}
#' \item{\code{cvup}}{Upper cross-validated errors.}
#' \item{\code{cvlo}}{Lower cross-validated errors.}
#' \item{\code{penalize.diagonal}}{Logical: whether the diagonal elements of \eqn{\mathbf{\Theta}} were penalized.}
#' \item{\code{diag.penalty.factor}}{The additional penalty multiplication factor for the diagonal elements of \eqn{\mathbf{\Theta}} when \code{'penalize.diagonal'} was returned as \code{'TRUE'}.}
#' @export
#' 
#' @author Yixiao Zeng \email{yixiao.zeng@@mail.mcgill.ca}, Celia M.T. Greenwood and Archer Yi Yang.
#' 
#' @examples
#' ## Simulate a dataset with response values missing completely at random (MCAR), 
#' ## the overall missing rate is around 10%.
#' sim.dat <- generateData(n = 300, p = 50, q = 20, rho = 0.1, missing.type = "MCAR")
#' tr <- 1:240  # training set indices
#' tst <- 241:300  # test set indices
#' X.tr <- sim.dat$X[tr, ]  # predictor matrix
#' Y.tr <- sim.dat$Z[tr, ]  # corrupted response matrix
#' 
#' \donttest{
#' ## Perform a five-fold cross-validation WITH specified 'lambda.Beta' and 'lambda.Theta'.
#' ## 'standardize' and 'standardize.response' are 'TRUE' by default.
#' lamB.vec <- 10^(seq(from = 0, to = -1, length.out = 5))
#' lamTht.vec <- 10^(seq(from = 0, to = -1, length.out = 5))
#' cvfit <- cv.missoNet(X = X.tr, Y = Y.tr, kfold = 5,
#'                      lambda.Beta = lamB.vec, lambda.Theta = lamTht.vec)
#' 
#' 
#' ## Perform a five-fold cross-validation WITHOUT specified 'lambda.Beta' and 'lambda.Theta'.
#' ## In this case, a grid of 'lambda.Beta' and 'lambda.Theta' values in a (hopefully) reasonable 
#' ## range will be computed and used by the program.
#' ## 
#' ## < This simplest command should be a good start for most analyses. >
#' cvfit <- cv.missoNet(X = X.tr, Y = Y.tr, kfold = 5)
#' 
#' 
#' ## Alternatively, compute the cross-validation folds in parallel on a cluster with 2 cores.
#' ## 
#' ## 'fit.1se = TRUE' tells the program to make additional estimations of the parameters at the 
#' ## largest value of 'lambda.Beta' / 'lambda.Theta' that gives the most regularized model such 
#' ## that the cross-validated error is within one standard error of the minimum.
#' cl <- parallel::makeCluster(min(parallel::detectCores()-1, 2))
#' cvfit <- cv.missoNet(X = X.tr, Y = Y.tr, kfold = 5, fit.1se = TRUE,
#'                      parallel = TRUE, cl = cl)
#' parallel::stopCluster(cl)
#' 
#' 
#' ## Use PRE-STANDARDIZED training data if you wish to compare the results with other softwares. 
#' ## There is no need for centering of variables.
#' X.tr.std <- scale(X.tr, center = FALSE, scale = apply(X.tr, 2, sd, na.rm = TRUE))
#' Y.tr.std <- scale(Y.tr, center = FALSE, scale = apply(Y.tr, 2, sd, na.rm = TRUE))
#' cvfit.std <- cv.missoNet(X = X.tr.std, Y = Y.tr.std, kfold = 5,
#'                          standardize = FALSE, standardize.response = FALSE)
#' 
#' 
#' ## Plot the (standardized) mean cross-validated errors in a heatmap.
#' plot(cvfit, type = "cv.heatmap")
#' 
#' ## Plot the (standardized) mean cross-validated errors in a 3D scatterplot.
#' plot(cvfit, type = "cv.scatter", plt.surf = TRUE)
#' 
#' 
#' ## Extract the estimates at "lambda.min".
#' Beta.hat1 <- cvfit$est.min$Beta
#' Theta.hat1 <- cvfit$est.min$Theta
#' 
#' ## Extract the estimates at "lambda.1se.Beta" (if 'fit.1se' = TRUE).
#' Beta.hat2 <- cvfit$est.1se.B$Beta
#' Theta.hat2 <- cvfit$est.1se.B$Theta
#' 
#' ## Extract the estimates at "lambda.1se.Theta" (if 'fit.1se' = TRUE).
#' Beta.hat3 <- cvfit$est.1se.Tht$Beta
#' Theta.hat3 <- cvfit$est.1se.Tht$Theta
#' 
#' 
#' ## Make predictions of response values on the test set.
#' newy1 <- predict(cvfit, newx = sim.dat$X[tst, ], s = "lambda.min")
#' newy2 <- predict(cvfit, newx = sim.dat$X[tst, ], s = "lambda.1se.Beta")  # 'fit.1se' = TRUE
#' newy3 <- predict(cvfit, newx = sim.dat$X[tst, ], s = "lambda.1se.Theta")  # 'fit.1se' = TRUE
#' }

cv.missoNet <- function(X, Y, kfold = 5, rho = NULL,
                        lambda.Beta = NULL, lambda.Theta = NULL,
                        lamBeta.min.ratio = NULL, lamTheta.min.ratio = NULL,
                        n.lamBeta = NULL, n.lamTheta = NULL,
                        lamBeta.scale.factor = 1, lamTheta.scale.factor = 1,
                        Beta.maxit = 1000, Beta.thr = 1e-04, eta = 0.8,
                        Theta.maxit = 1000, Theta.thr = 1e-04, eps = 1e-08,
                        penalize.diagonal = NULL, diag.penalty.factor = NULL,
                        standardize = TRUE, standardize.response = TRUE,
                        fit.1se = FALSE, fit.relax = FALSE,
                        permute = TRUE, with.seed = NULL,
                        parallel = FALSE, cl = NULL, verbose = 1) {
  if (verbose > 0) { cat("\n======================= cv.missoNet =======================\n
- Parameter initialization ...\n\n") }
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  if (permute) {
    set.seed(with.seed)
    ind <- sample(n, replace = FALSE)
  } else { ind <- 1:n }
  foldid <- unlist(lapply(1:kfold, function(x) { rep(x, length((1 + floor((x - 1) * n/kfold)):floor(x * n/kfold))) }))
  names(ind) <- paste0("fold.", foldid)
  
  init.obj <- InitParams(X = X[ind, ], Y = Y[ind, ], rho = rho, kfold = kfold, foldid = foldid,
                         Theta.maxit = Theta.maxit, Theta.thr = Theta.thr, eps = eps, 
                         penalize.diagonal = penalize.diagonal, diag.pf = diag.penalty.factor,
                         standardize = standardize, standardize.response = standardize.response)
  
  lambda.obj <- InitLambda(lamB = lambda.Beta, lamTh = lambda.Theta, n.tr = floor(n * (kfold - 1)/kfold),
                           init.obj = init.obj, n.lamB = n.lamBeta, n.lamTh = n.lamTheta,
                           lamB.min.ratio = lamBeta.min.ratio, lamTh.min.ratio = lamTheta.min.ratio,
                           lamB.scale.factor = lamBeta.scale.factor, lamTh.scale.factor = lamTheta.scale.factor)
  lamTh.vec <- lambda.obj$lamTh.vec
  lamB.vec <- lambda.obj$lamB.vec
  
  ################################################################################
  # Cross-validation
  ################################################################################
  if (verbose > 0) { cat("--------------------- Cross-validation --------------------\n\n") }
  if (!parallel) {
    err <- matrix(0, kfold, length(lamTh.vec))
    for (k in 1:kfold) {
      if (verbose > 0) { cat(sprintf("Fold: %d/%d\n", k, kfold)) }
      
      foldind <- ind[(1 + floor((k - 1) * n/kfold)):floor(k * n/kfold)]
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
      rho.mat.1 <- t(matrix(rep(1 - rho.vec, p), q, p))  ## pxq
      rho.mat.2 <- matrix(1 - rho.vec, q, 1) %*% matrix(1 - rho.vec, 1, q)
      diag(rho.mat.2) <- 1 - rho.vec  ## qxq
      
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
      info.update$B.init <- init.obj$B.init * init.obj$sdx   ## initialize B on the standardized scale
      Beta.thr.rescale <- Beta.thr * sum(abs(info.update$B.init))
      E.tr <- Y.tr - X.tr %*% info.update$B.init
      info.update$residual.cov <- getResidual(E = E.tr, n = n.tr, rho.mat = rho.mat.2, eps = eps)
      
      if (verbose == 1) { pb <- txtProgressBar(min = 0, max = length(lamTh.vec), style = 3, width = 50, char = "=") }
      for (i in 1:length(lamTh.vec)) {
        info.update$B.init <- update.missoNet(lamTh = lamTh.vec[i], lamB = lamB.vec[i],
                                              Beta.maxit = Beta.maxit, Beta.thr = Beta.thr.rescale,
                                              Theta.maxit = Theta.maxit, Theta.thr = Theta.thr,
                                              verbose = verbose, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                                              info = info, info.update = info.update, under.cv = TRUE)
        Beta.thr.rescale <- Beta.thr * sum(abs(info.update$B.init))
        E.tr <- Y.tr - X.tr %*% info.update$B.init
        info.update$residual.cov <- getResidual(E = E.tr, n = n.tr, rho.mat = rho.mat.2, eps = eps)
        
        E.va.sq <- (Y.va - X.va %*% info.update$B.init)^2
        err[k, i] <- mean(E.va.sq, na.rm = TRUE)
        if (verbose == 1) { setTxtProgressBar(pb, i) }
      }
      if (verbose == 1) { close(pb) }
    }

  } else {
    if (verbose > 0) {
      cat("- Parallel execution on", length(cl), "CPU cores ...\n\n")
      pbapply::pboptions(type = "txt", style = 3, char = "=", txt.width = 50, use_lb = TRUE, nout = kfold)
    } else {pbapply::pboptions(type = "none", use_lb = TRUE)}
   
    par.out <- pbapply::pblapply(1:kfold, function(k) {
      parWrapper(k = k, X = X, Y = Y, init.obj = init.obj, rho = rho, ind = ind, kfold = kfold, lamTh.vec = lamTh.vec, lamB.vec = lamB.vec,
                 Beta.maxit = Beta.maxit, Beta.thr = Beta.thr, Theta.maxit = Theta.maxit, Theta.thr = Theta.thr, eps = eps, eta = eta)
    }, cl = cl)
    
    err <- do.call("rbind", par.out)
  }
  if (verbose > 0) { cat("\n-----------------------------------------------------------\n\n") }
  
  err.cv <- colSums(err)/kfold
  err.sd <- apply(err, 2, sd)/sqrt(kfold)
  err.up <- err.cv + err.sd
  err.low <- err.cv - err.sd
  
  cv.min <- which.min(err.cv)
  lamTh.min <- lamTh.vec[cv.min]
  lamB.min <- lamB.vec[cv.min]
  boundaryCheck(lambda.Theta = lambda.Theta, lambda.Beta = lambda.Beta, 
                lamTh.vec = lamTh.vec, lamB.vec = lamB.vec, 
                lamTh.min = lamTh.min, lamB.min = lamB.min, margin = 0.1)
  
  if (verbose > 0) { cat("- Fittig with `lambda.min` ...\n\n") }
  out.min <- update.missoNet(X = X, Y = Y, lamTh = lamTh.min, lamB = lamB.min,
                             Beta.maxit = Beta.maxit * 10, Beta.thr = Beta.thr * 0.01,
                             Theta.maxit = Theta.maxit * 10, Theta.thr = Theta.thr * 0.01,
                             verbose = verbose, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                             info = NULL, info.update = NULL, init.obj = init.obj, under.cv = FALSE)
  out.min$Beta <- sweep(out.min$Beta/init.obj$sdx, 2, init.obj$sdy, `*`)    ## convert back to the original scale
  out.min$mu <- as.numeric(init.obj$my - crossprod(out.min$Beta, init.obj$mx))
  out.min$lambda.Beta <- lamB.min
  out.min$lambda.Theta <- lamTh.min
  relax.net <- NULL
  if (fit.relax) {
    if (verbose > 0) { cat("  -- Relaxed network\n\n") }
    relax.net <- relax.glasso(X = X, Y = Y, init.obj = init.obj, est = out.min, eps = eps,
                              Theta.thr = Theta.thr * 0.01, Theta.maxit = Theta.maxit * 10)
  }
  out.min$relax.net <- relax.net
  
  outB.1se <- NULL
  outTh.1se <- NULL
  if (fit.1se) {
    if (verbose > 0) { 
      cat("------------------------------\n
- Fiting with `lambda.1se` ...\n
  -- `lambda.1se.Beta`\n\n") }
    new.lamB.vec <- lamB.vec[lamTh.vec == lamTh.min]
    new.err.cv <- err.cv[lamTh.vec == lamTh.min]
    lamB.1se <- max(new.lamB.vec[new.err.cv <= err.up[cv.min]])
    
    if (lamB.1se != lamB.min) {
      outB.1se <- update.missoNet(X = X, Y = Y, lamTh = lamTh.min, lamB = lamB.1se,
                                  Beta.maxit = Beta.maxit * 10, Beta.thr = Beta.thr * 0.01,
                                  Theta.maxit = Theta.maxit * 10, Theta.thr = Theta.thr * 0.01,
                                  verbose = verbose, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                                  info = NULL, info.update = NULL, init.obj = init.obj, under.cv = FALSE)
      outB.1se$Beta <- sweep(outB.1se$Beta/init.obj$sdx, 2, init.obj$sdy, `*`)
      outB.1se$mu <- as.numeric(init.obj$my - crossprod(outB.1se$Beta, init.obj$mx))
      outB.1se$lambda.Beta <- lamB.1se
      outB.1se$lambda.Theta <- lamTh.min
      relax.net <- NULL
      if (fit.relax) {
        if (verbose > 0) { cat("    --- Relaxed network\n\n") }
        relax.net <- relax.glasso(X = X, Y = Y, init.obj = init.obj, est = outB.1se, eps = eps,
                                  Theta.thr = Theta.thr * 0.01, Theta.maxit = Theta.maxit * 10)
      }
      outB.1se$relax.net <- relax.net
    } else {
      warning("\n`lambda.1se.Beta` = `lambda.min.Beta`, please provide a finer grid for `lambda.Beta` by increasing the number of values.\n")
    }
    
    if (verbose > 0) { cat("  -- `lambda.1se.Theta`\n\n") }
    new.lamTh.vec <- lamTh.vec[lamB.vec == lamB.min]
    new.err.cv <- err.cv[lamB.vec == lamB.min]
    lamTh.1se <- max(new.lamTh.vec[new.err.cv <= err.up[cv.min]])
    
    if (lamTh.1se != lamTh.min) {
      outTh.1se <- update.missoNet(X = X, Y = Y, lamTh = lamTh.1se, lamB = lamB.min,
                                  Beta.maxit = Beta.maxit * 10, Beta.thr = Beta.thr * 0.01,
                                  Theta.maxit = Theta.maxit * 10, Theta.thr = Theta.thr * 0.01,
                                  verbose = verbose, eps = eps, eta = eta, diag.pf = init.obj$diag.pf,
                                  info = NULL, info.update = NULL, init.obj = init.obj, under.cv = FALSE)
      outTh.1se$Beta <- sweep(outTh.1se$Beta/init.obj$sdx, 2, init.obj$sdy, `*`)
      outTh.1se$mu <- as.numeric(init.obj$my - crossprod(outTh.1se$Beta, init.obj$mx))
      outTh.1se$lambda.Beta <- lamB.min
      outTh.1se$lambda.Theta <- lamTh.1se
      relax.net <- NULL
      if (fit.relax) {
        if (verbose > 0) { cat("    --- Relaxed network\n\n") }
        relax.net <- relax.glasso(X = X, Y = Y, init.obj = init.obj, est = outTh.1se, eps = eps,
                                  Theta.thr = Theta.thr * 0.01, Theta.maxit = Theta.maxit * 10)
      }
      outTh.1se$relax.net <- relax.net
    } else {
      warning("\n`lambda.1se.Theta` = `lambda.min.Theta`, please provide a finer grid for `lambda.Theta` by increasing the number of values.\n")
    }
  }
  
  if (verbose > 0) { cat("========================= FINISHED ========================\n\n") }
  
  cv <- list(est.min = out.min, est.1se.B = outB.1se, est.1se.Tht = outTh.1se, rho = init.obj$rho.vec, fold.index = ind,
             lambda.Beta.vec = lamB.vec, lambda.Theta.vec = lamTh.vec, cvm = err.cv, cvup = err.up, cvlo = err.low,
             penalize.diagonal = init.obj$penalize.diagonal, diag.penalty.factor = init.obj$diag.pf)
  class(cv) <- c("cv.missoNet", class(cv))
  return(cv)
}

