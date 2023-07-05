<!-- badges: start -->
[![R-CMD-check](https://github.com/yixiao-zeng/missoNet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yixiao-zeng/missoNet/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/missoNet)](https://CRAN.R-project.org/package=missoNet)
<!-- badges: end -->

# missoNet: Missingness in Multi-Task Regression with Network Estimation

`missoNet` is a novel approach to fitting penalized multi-task regression models, which are used to 
estimate the coefficients of predictor variables for multiple correlated tasks/response variables. 
The package achieves this by simultaneously estimating the regression coefficients and the conditional 
response network structure given all predictors, using penalized maximum likelihood in an undirected 
conditional Gaussian graphical model. In contrast to most penalized multi-task regression methods, such 
as conditional graphical lasso, `missoNet` is capable of obtaining estimates even when the response data 
is corrupted by missing values. The method is based on convex optimization, which provides both theoretical 
and computational benefits, and returns solutions that are comparable to the estimates obtained without 
any missing values.

The package provides an integrated set of core routines including 1) generation of simulation data; 2) model fitting and 
cross-validation; 3) visualization of results; 4) predictions in new data. The function arguments are specified 
in the same style as those of `glmnet`, making it easy for experienced users to get started.


## Installation

To install the package `missoNet` from CRAN, type the following command in the R console:

```{r}
install.packages("missoNet")
```

Or install the development version of `missoNet` from GitHub:

```r
if(!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("yixiao-zeng/missoNet", build_vignettes = TRUE)
```


## An example for getting started

An example of how to use the package:

```r
# Simulate a dataset with response values missing completely at random (MCAR), 
# the overall missing rate is around 10%.
sim.dat <- generateData(n = 300, p = 50, q = 20, rho = 0.1, missing.type = "MCAR")
tr <- 1:240  # training set indices
tst <- 241:300  # test set indices
X.tr <- sim.dat$X[tr, ]  # predictor matrix
Y.tr <- sim.dat$Z[tr, ]  # corrupted response matrix

# Perform a five-fold cross-validation on the training data.
cvfit <- cv.missoNet(X = X.tr, Y = Y.tr, kfold = 5)

# Alternatively, compute the cross-validation folds in parallel.
cl <- parallel::makeCluster(min(parallel::detectCores()-1, 3))
cvfit <- cv.missoNet(X = X.tr, Y = Y.tr, kfold = 5,
                     parallel = TRUE, cl = cl)
parallel::stopCluster(cl)

# Plot the standardized mean cross-validated errors in a heatmap.
plot(cvfit)

# Extract the estimates at "lambda.min" that gives the minimum cross-validated error.
Beta_hat <- cvfit$est.min$Beta
Theta_hat <- cvfit$est.min$Theta

# Make predictions of response values on the test set.
newy <- predict(cvfit, newx = sim.dat$X[tst, ], s = "lambda.min")
```


## Learn more

See the vignette for more detailed information.

```r
vignette("missoNet")
```


## References

- [Cocolasso for high dimensional error-in-variables regression](https://arxiv.org/pdf/1510.07123.pdf)

- [High-dimensional regression with noisy and missing data: Provable guarantees with nonconvexity](https://arxiv.org/pdf/1109.3714.pdf)

