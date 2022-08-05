<!-- badges: start -->
[![R-CMD-check](https://github.com/yixiao-zeng/missoNet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yixiao-zeng/missoNet/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# missoNet: multi-task regression and conditional network estimation with missing data.

The package **missoNet** enables users to fit multi-task Gaussian regression models 
to estimate covariate effects by leveraging the relatedness information of multiple responses (outputs). 
Meanwhile, the algorithm can learn the conditional network structure among outputs in an undirected Gaussian graphical model.

Different from existing techniques for multivariate regression/multi-task learning, **missoNet** allows missing data in the 
output matrix, it enjoys the theoretical and computational benefits of convexity and returns solutions that are 
comparable/close to the clean data estimates.

The package includes methods for cross-validation, and functions for prediction and plotting. It has function arguments in the 
same style as those of **glmnet**, making it easy for experienced users to get started with.


## Installation

Install the development version of `missoNet` from GitHub:

```r
if(!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("yixiao-zeng/missoNet")
```


## An example for getting started

An example of how to use the package:

```r
## generate a simulated dataset with overall 
## missing rate = 0.1, missing mechanism = "MCAR"
sim.dat <- generateData(n = 300, p = 50, q = 20, rho = 0.1, missing.type = "MCAR")
tr <- 1:240  ## training set
va <- 241:300  ## validation set

## use `cv.missoNet` to do a five-fold cross-validation 
cvfit <- cv.missoNet(X = sim.dat$X[tr, ], Y = sim.dat$Z[tr, ], kfold = 5)

## or train the model in parallel
library(snowfall)
cvfit <- cv.missoNet(X = sim.dat$X[tr, ], Y = sim.dat$Z[tr, ], kfold = 5,
                     parallel = TRUE, cpus = min(parallel::detectCores()-1, 5)) 

## plot the heatmap of the cross-validation error
plot(cvfit)

## parameters estimated `lambda.min` that gives the smallest CV error
B_hat <- cvfit$est.min$Beta
Tht_hat <- cvfit$est.min$Theta

## make prediction
newy <- predict(cvfit, newx = sim.dat$X[va, ], s = "lambda.min")
```


## Learn more

See the vignette for more detailed information.

```r
vignette("missoNet")
```


## References

- [Cocolasso for high dimensional error-in-variables regression](https://arxiv.org/pdf/1510.07123.pdf)

- [High-dimensional regression with noisy and missing data: Provable guarantees with nonconvexity](https://arxiv.org/pdf/1109.3714.pdf)

