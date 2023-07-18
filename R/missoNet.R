#' Fit a series of missoNet models with user-supplied regularization parameters for the lasso penalties
#'
#' This function fits the conditional graphical lasso models to datasets with missing response values. 
#' \sQuote{\code{missoNet}} computes the regularization path for the lasso penalties sequentially along the
#' bivariate regularization parameter sequence \eqn{\{(\lambda_B, \lambda_\Theta)\}} provided by the user.
#' 
#' \sQuote{\code{missoNet}} is the main model-fitting function which is specifically proposed to fit the conditional 
#' graphical lasso models / penalized multi-task Gaussian regressions to (corrupted) datasets with response values missing at random (MAR).
#' To facilitate the interpretation of the model, let's temporarily assume that there are no missing values 
#' in the data used to fit the model. Suppose we have \eqn{n} observations of both a \eqn{p}-variate predictor \eqn{X \in \mathcal{R}^p}
#' and a \eqn{q}-variate response \eqn{Y \in \mathcal{R}^q}, for the \eqn{i}th sample (\eqn{i = 1,...,n}), 
#' \sQuote{\code{missoNet}} assumes the model
#' \deqn{Y_i = \mu + X_i\mathbf{B} + E_i,\ \ E_i \sim \mathcal{MVN}(0_q, (\mathbf{\Theta})^{-1}),}
#' where \eqn{Y_i \in \mathcal{R}^{1\times q}} and \eqn{X_i \in \mathcal{R}^{1\times p}} are one 
#' realization of the \eqn{q} responses and the \eqn{p} predictors, respectively. 
#' \eqn{E_i \in \mathcal{R}^{1\times q}} is an error vector drawn from a multivariate Gaussian distribution. 
#' 
#' The regression coefficient matrix \eqn{\mathbf{B} \in \mathcal{R}^{p\times q}} that mapping predictors to responses and 
#' the precision (inverse covariance) matrix \eqn{\mathbf{\Theta} \in \mathcal{R}^{q\times q}} that revealing the 
#' responses' conditional dependencies are the parameters to be estimated by solving a penalized MLE problem
#' \deqn{(\hat{\mathbf{\Theta}},\hat{\mathbf{B}}) = {\mathrm{argmin}}_{\mathbf{\Theta} \succeq 0,\ \mathbf{B}}\ 
#' g(\mathbf{\Theta},\mathbf{B}) + \lambda_{\Theta}(\|\mathbf{\Theta}\|_{1,\mathrm{off}} + 1_{n\leq \mathrm{max}(p,q)} \|\mathbf{\Theta}\|_{1,\mathrm{diag}}) + \lambda_{B}\|\mathbf{B}\|_1,}
#' where 
#' \deqn{g(\mathbf{\Theta},\mathbf{B}) = \mathrm{tr}\left[\frac{1}{n}(\mathbf{Y} - \mathbf{XB})^\top(\mathbf{Y} - \mathbf{XB}) \mathbf{\Theta}\right] 
#' - \mathrm{log}|\mathbf{\Theta}|.}
#' The response matrix \eqn{\mathbf{Y} \in \mathcal{R}^{n\times q}} has \eqn{i}th row (\eqn{Y_i - \frac{1}{n}\sum_{j=1}^n Y_j}), 
#' and the predictor matrix \eqn{\mathbf{X} \in \mathcal{R}^{n\times p}} has \eqn{i}th row (\eqn{X_i - \frac{1}{n}\sum_{j=1}^n X_j}). 
#' The intercept \eqn{\mu \in \mathcal{R}^{1\times q}} is canceled out because of centering of the data matrices \eqn{\mathbf{Y}} and \eqn{\mathbf{X}}. 
#' \eqn{1_{n\leq \mathrm{max}(p,q)}} denotes the indicator function for whether penalizing the diagonal elements of \eqn{\mathbf{\Theta}} or not. 
#' When \eqn{n\leq \mathrm{max}(p,q)}, a global minimizer of the objective function defined above does not exist without the diagonal penalization.
#' 
#' Missingness in real data is inevitable. In this instance, the estimates based only on complete cases are likely to be biased, 
#' and the objective function is likely to no longer be a biconvex optimization problem. In addition, many algorithms cannot be directly employed since they 
#' require complete datasets as inputs. \sQuote{\code{missoNet}} aims to handle the specific situation where the response matrix \eqn{\mathbf{Y}} contains values that 
#' are missing at random (MAR. Please refer to the vignette or other resources for more information about the differences between MAR, missing completely at 
#' random (MCAR) and missing not at random (MNAR)). As it should be, \sQuote{\code{missoNet}} is also applicable to datasets with MCAR response values or without any missing values. 
#' The method provides a unified framework for automatically solving a convex modification of the multi-task learning problem defined above, 
#' using corrupted datasets. Moreover, \sQuote{\code{missoNet}} enjoys the theoretical and computational benefits of convexity and returns 
#' solutions that are comparable/close to the clean conditional graphical lasso estimates. Please refer to the original manuscript (coming soon) for more details of our method.
#' 
#' @param X Numeric predictor matrix (\eqn{n\times p}): columns correspond to predictor variables and rows correspond to samples. Missing values are not allowed. There is no need for centering or scaling of the variables. \code{'X'} should not include a column of ones for an intercept.
#' @param Y Numeric response matrix (\eqn{n\times q}): columns correspond to response variables and rows correspond to samples. Missing values should be coded as either \code{'NA'}s or \code{'NaN'}s. There is no need for centering or scaling of the variables.
#' @param lambda.Beta A scalar or a numeric vector: a user-supplied sequence of non-negative value(s) for \{\eqn{\lambda_B}\} used to penalize the elements of the coefficient matrix \eqn{\mathbf{B}}. Note that the values will be sequentially visited in the given orders as inputs to the regularization parameter sequence \eqn{\{(\lambda_B, \lambda_\Theta)\}}; \code{'lambda.Beta'} must have the same length as \code{'lambda.Theta'}.
#' @param lambda.Theta A scalar or a numeric vector: a user-supplied sequence of non-negative value(s) for \{\eqn{\lambda_\Theta}\} used to penalize the (off-diagonal) elements of the precision matrix \eqn{\mathbf{\Theta}}. Note that the values will be sequentially visited in the given orders as inputs to the regularization parameter sequence \eqn{\{(\lambda_B, \lambda_\Theta)\}}; \code{'lambda.Theta'} must have the same length as \code{'lambda.Beta'}.
#' @param rho (Optional) A scalar or a numeric vector of length \eqn{q}: the elements are user-supplied probabilities of missingness for the response variables. The default is \code{'rho = NULL'} and the program will compute the empirical missing rates for each of the columns of \code{'Y'} and use them as the working missing probabilities. The default setting should suffice in most cases; misspecified missing probabilities would introduce biases into the model.
#' @param Beta.maxit The maximum number of iterations of the fast iterative shrinkage-thresholding algorithm (FISTA) for updating \eqn{\hat{\mathbf{B}}}. The default is \code{'Beta.maxit = 10000'}.
#' @param Beta.thr The convergence threshold of the FISTA algorithm for updating \eqn{\hat{\mathbf{B}}}; the default is \code{'Beta.thr = 1.0E-8'}. Iterations stop when the absolute parameter change is less than (\code{'Beta.thr'} \code{*} \code{sum(abs(}\eqn{\hat{\mathbf{B}}}\code{))}).
#' @param eta The backtracking line search shrinkage factor; the default is \code{'eta = 0.8'}. Most users will be able to use the default value, some experienced users may want to tune \code{'eta'} according to the properties of a specific dataset for a faster convergence of the FISTA algorithm. Note that \code{'eta'} must be in (0, 1).
#' @param Theta.maxit The maximum number of iterations of the \sQuote{\code{\link{glasso}}} algorithm for updating \eqn{\hat{\mathbf{\Theta}}}. The default is \code{'Theta.maxit = 10000'}.
#' @param Theta.thr The convergence threshold of the \sQuote{\code{\link{glasso}}} algorithm for updating \eqn{\hat{\mathbf{\Theta}}}; the default is \code{'Theta.thr = 1.0E-8'}. Iterations stop when the average absolute parameter change is less than (\code{'Theta.thr'} \code{*} \code{ave(abs(offdiag(}\eqn{\hat{\mathbf{\Sigma}}}\code{)))}), where \eqn{\hat{\mathbf{\Sigma}}} denotes the empirical working covariance matrix.
#' @param eps A numeric tolerance level for the L1 projection of the empirical covariance matrix; the default is \code{'eps = 1.0E-8'}. The empirical covariance matrix will be projected onto a L1 ball to have \code{min(eigen(}\eqn{\hat{\mathbf{\Sigma}}}\code{)$value)} == \code{'eps'}, if any of the eigenvalues is less than the specified tolerance. Most users will be able to use the default value.
#' @param penalize.diagonal Logical: should the diagonal elements of \eqn{\mathbf{\Theta}} be penalized? The default is \code{'TRUE'}.
#' @param diag.penalty.factor Numeric: a separate penalty multiplication factor for the diagonal elements of \eqn{\mathbf{\Theta}} when \code{'penalize.diagonal = TRUE'}. \eqn{\lambda_\Theta} is multiplied by this number to allow a differential shrinkage of the diagonal elements. The default is \code{'NULL'} and the program will guess a value based on an initial estimate of \eqn{\mathbf{\Theta}}. This factor could be \code{'0'} for no shrinkage (equivalent to \code{'penalize.diagonal = FALSE'}) or \code{'1'} for an equal shrinkage.
#' @param standardize Logical: should the columns of \code{'X'} be standardized so each has unit variance? The default is \code{'TRUE'}. The estimated results will always be returned on the original scale. If \code{'X'} has been standardized prior to fitting the model, you might not wish to standardize it inside the algorithm.
#' @param standardize.response Logical: should the columns of \code{'Y'} be standardized so each has unit variance? The default is \code{'TRUE'}. The estimated results will always be returned on the original scale. If \code{'Y'} has been standardized prior to fitting the model, you might not wish to standardize it inside the algorithm.
#' @param fit.relax Logical: the default is \code{'FALSE'}. If \code{'TRUE'}, the program will re-estimate the edges in the active set (i.e. nonzero off-diagonal elements) of the network structure \eqn{\hat{\mathbf{\Theta}}} without penalization (\eqn{\lambda_\Theta=0}). This debiased estimate of \eqn{\mathbf{\Theta}} could be useful for some interdependency analyses. WARNING: there may be convergence issues if the empirical covariance matrix is not of full rank (e.g. \eqn{n < q)}).
#' @param parallel Logical: the default is \code{'FALSE'}. If \code{'TRUE'}, the program uses clusters to fit models with each element of the \eqn{\lambda} sequence \eqn{\{(\lambda_B, \lambda_\Theta)\}} in parallel. Must register parallel clusters beforehand, see examples below.
#' @param cl A cluster object created by \sQuote{\code{parallel::makeCluster}} for parallel evaluations. This is only needed when \code{'parallel = TRUE'}.
#' @param verbose Value of \code{'0'}, \code{'1'} or \code{'2'}. \code{'verbose = 0'} -- silent; \code{'verbose = 1'} (the default) -- limited tracing with progress bars; \code{'verbose = 2'} -- detailed tracing. Note that displaying the progress bars slightly increases the computation overhead compared to the silent mode. The detailed tracing should be used for monitoring progress only when the program runs extremely slowly, and it is not supported under \code{'parallel = TRUE'}.
#'
#' @return This function returns a \code{'list'} consisting of the following components:
#' \item{\code{est.list}}{A named \code{'list'} storing the lists of results estimated at each of the \eqn{\lambda} pairs, (\eqn{\lambda_B}, \eqn{\lambda_\Theta}). Each sub-\code{'list'} contains:
#'   \itemize{
#'     \item \code{Beta}: the penalized estimate of the regression coefficient matrix \eqn{\hat{\mathbf{B}}} (\eqn{p\times q}).
#'     \item \code{Theta}: the penalized estimate of the precision matrix \eqn{\hat{\mathbf{\Theta}}} (\eqn{q\times q}).
#'     \item \code{mu}: a vector of length \eqn{q} storing the model intercept \eqn{\hat{\mu}}.
#'     \item \code{lambda.Beta}: the value of \eqn{\lambda_B} used to fit the model.
#'     \item \code{lambda.Theta}: the value of \eqn{\lambda_\Theta} used to fit the model.
#'     \item \code{relax.net}: the relaxed (debiased) estimate of the conditional network structure \eqn{\hat{\mathbf{\Theta}}_\mathrm{rlx}} (\eqn{q\times q}) if \code{'fit.relax = TRUE'} when calling \sQuote{\code{missoNet}}.
#'   }
#' }
#' \item{\code{rho}}{A vector of length \eqn{q} storing the working missing probabilities for the \eqn{q} response variables.}
#' \item{\code{penalize.diagonal}}{Logical: whether the diagonal elements of \eqn{\mathbf{\Theta}} were penalized.}
#' \item{\code{diag.penalty.factor}}{The additional penalty multiplication factor for the diagonal elements of \eqn{\mathbf{\Theta}} when \code{'penalize.diagonal'} was returned as \code{'TRUE'}.}
#' @export
#'
#' @author Yixiao Zeng \email{yixiao.zeng@@mail.mcgill.ca}, Celia M.T. Greenwood and Archer Yi Yang.
#'
#' @examples
#' ## Simulate a dataset with response values missing completely at random (MCAR), 
#' ## the overall missing rate is around 10%.
#' set.seed(123)  # reproducibility
#' sim.dat <- generateData(n = 300, p = 50, q = 20, rho = 0.1, missing.type = "MCAR")
#' tr <- 1:240  # training set indices
#' tst <- 241:300  # test set indices
#' X.tr <- sim.dat$X[tr, ]  # predictor matrix
#' Y.tr <- sim.dat$Z[tr, ]  # corrupted response matrix
#' 
#' \donttest{
#' ## Fit one missoNet model with two scalars for 'lambda.Beta' and 'lambda.Theta'.
#' fit1 <- missoNet(X = X.tr, Y = Y.tr, lambda.Beta = 0.1, lambda.Theta = 0.2)
#' 
#' 
#' ## Fit a series of missoNet models with the lambda pairs := (lambda.Beta, lambda.Theta)
#' ## sequentially extracted from the 'lambda.Beta' and 'lambda.Theta' vectors, note that the 
#' ## two vectors must have the same length.
#' lamB.vec <- 10^(seq(from = 0, to = -1, length.out = 5))
#' lamTht.vec <- rep(0.1, 5)
#' fit2 <- missoNet(X = X.tr, Y = Y.tr, lambda.Beta = lamB.vec, lambda.Theta = lamTht.vec)
#' 
#' 
#' ## Parallelization on a cluster with two cores.
#' cl <- parallel::makeCluster(2)
#' fit2 <- missoNet(X = X.tr, Y = Y.tr, lambda.Beta = lamB.vec, lambda.Theta = lamTht.vec, 
#'                  parallel = TRUE, cl = cl)
#' parallel::stopCluster(cl)
#' 
#' 
#' ## Extract the estimates at ('lamB.vec[1]', 'lamTht.vec[1]').
#' ## The estimates at the subsequent lambda pairs could be accessed in the same way.
#' Beta.hat <- fit2$est.list[[1]]$Beta
#' Theta.hat <- fit2$est.list[[1]]$Theta
#' lambda.Beta <- fit2$est.list[[1]]$lambda.Beta  # equal to 'lamB.vec[1]'
#' lambda.Theta <- fit2$est.list[[1]]$lambda.Theta  # equal to 'lamTht.vec[1]'
#' 
#' 
#' ## Fit a series of missoNet models using PRE-STANDARDIZED training data
#' ## if you wish to compare the results with other softwares. 
#' X.tr.std <- scale(X.tr, center = TRUE, scale = TRUE)
#' Y.tr.std <- scale(Y.tr, center = TRUE, scale = TRUE)
#' fit3 <- missoNet(X = X.tr.std, Y = Y.tr.std, lambda.Beta = lamB.vec, lambda.Theta = lamTht.vec,
#'                  standardize = FALSE, standardize.response = FALSE)
#' }

missoNet <- function(X, Y, lambda.Beta, lambda.Theta, rho = NULL,
                     Beta.maxit = 1e4, Beta.thr = 1e-08, eta = 0.8,
                     Theta.maxit = 1e4, Theta.thr = 1e-08, eps = 1e-08,
                     penalize.diagonal = TRUE, diag.penalty.factor = NULL,
                     standardize = TRUE, standardize.response = TRUE,
                     fit.relax = FALSE, parallel = FALSE, cl = NULL, verbose = 1) {
  if (length(lambda.Beta) != length(lambda.Theta)) {
    stop("`lambda.Beta` and `lambda.Theta` should be equal in length.")
  }
  
  if (verbose > 0) { cat("\n========================= missoNet ========================\n
- Model initialization ...\n\n") }
  init.obj <- InitParams(X = X, Y = Y, rho = rho, under.cv = FALSE, lamB.vec = lambda.Beta,
                         eps = eps, penalize.diag = penalize.diagonal, diag.pf = diag.penalty.factor,
                         standardize = standardize, standardize.response = standardize.response)
  B.init.list <- init.obj$B.init
  init.obj$B.init <- NULL
  
  if (verbose > 0) { cat("-----------------------------------------------------------\n
- Fittig with the user-supplied `lambda` pair(s) ...\n\n") }
  if (!parallel) {
    fit_list <- vector("list", length(lambda.Theta))
    names(fit_list) <- paste0(1:length(lambda.Theta), ": lamB=", sprintf("%.3f", lambda.Beta), " lamTht=", sprintf("%.3f", lambda.Theta))
    if (verbose == 1) { pb <- txtProgressBar(min = 0, max = length(lambda.Theta), style = 3, width = 50, char = "=") }
    for (i in 1:length(lambda.Theta)) {
      fit_list[[i]] <- fitWrapper(X = X, Y = Y, lambda.Theta = lambda.Theta[i], lambda.Beta = lambda.Beta[i],
                                  Beta.maxit = Beta.maxit, Beta.thr = Beta.thr, eta = eta,
                                  Theta.maxit = Theta.maxit, Theta.thr = Theta.thr, eps = eps,
                                  penalize.diagonal = penalize.diagonal, verbose = verbose, fit.relax = fit.relax,
                                  init.obj = init.obj, B.init = B.init.list[[i]])
      if (verbose == 1) { setTxtProgressBar(pb, i) }
    }
    if (verbose == 1) { close(pb) }
  } else {
    if (verbose > 0) {
      cat("  - Parallel execution on", length(cl), "CPU cores ...\n\n")
      pbapply::pboptions(type = "txt", style = 3, char = "=", txt.width = 50, use_lb = TRUE, nout = min(length(lambda.Theta), 100))
    } else {pbapply::pboptions(type = "none", use_lb = TRUE)}
    
    fit_list <- pbapply::pblapply(1:length(lambda.Theta), function(i) {
      fitWrapper(X = X, Y = Y, lambda.Theta = lambda.Theta[i], lambda.Beta = lambda.Beta[i],
                 Beta.maxit = Beta.maxit, Beta.thr = Beta.thr, eta = eta,
                 Theta.maxit = Theta.maxit, Theta.thr = Theta.thr, eps = eps,
                 penalize.diagonal = penalize.diagonal, verbose = 0, fit.relax = fit.relax,
                 init.obj = init.obj, B.init = B.init.list[[i]])
    }, cl = cl)
    
    names(fit_list) <- paste0(1:length(lambda.Theta), ": lamB=", sprintf("%.3f", lambda.Beta), " lamTht=", sprintf("%.3f", lambda.Theta))
  }
  if (verbose > 0) { cat("\n========================= FINISHED ========================\n\n") }
  
  return(list(est.list = fit_list, rho = init.obj$rho.vec,
              penalize.diagonal = penalize.diagonal, diag.penalty.factor = init.obj$diag.pf))
}


fitWrapper <- function(X, Y, lambda.Theta, lambda.Beta,
                       Beta.maxit, Beta.thr, eta,
                       Theta.maxit, Theta.thr, eps,
                       penalize.diagonal, verbose, fit.relax,
                       init.obj, B.init) {
  
  fit <- update.missoNet(X = X, Y = Y, lamTh = lambda.Theta, lamB = lambda.Beta,
                         Beta.maxit = Beta.maxit, Beta.thr = Beta.thr,
                         Theta.maxit = Theta.maxit, Theta.thr = Theta.thr,
                         verbose = verbose, eps = eps, eta = eta,
                         penalize.diag = penalize.diagonal, diag.pf = init.obj$diag.pf,
                         info = NULL, info.update = NULL, under.cv = FALSE,
                         init.obj = init.obj, B.init = B.init)
  fit$Beta <- sweep(fit$Beta/init.obj$sdx, 2, init.obj$sdy, `*`)    ## convert back to the original scale
  fit$mu <- as.numeric(init.obj$my - crossprod(fit$Beta, init.obj$mx))
  fit$lambda.Beta <- lambda.Beta
  fit$lambda.Theta <- lambda.Theta
  relax.net <- NULL
  if (fit.relax) {
    relax.net <- relax.glasso(X = X, Y = Y, init.obj = init.obj, est = fit, eps = eps,
                              Theta.thr = Theta.thr, Theta.maxit = Theta.maxit)
  }
  fit$relax.net <- relax.net
  
  return(fit)
}

