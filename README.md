<!-- badges: start -->
[![R-CMD-check](https://github.com/yixiao-zeng/missoNet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yixiao-zeng/missoNet/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# missoNet: multi-task regression and conditional network estimation with missing data.

The package **missoNet** enables users to fit multi-task gaussian regression models 
and estimate covariate effects by leveraging the relatedness information of multiple responses (outputs). 
Meanwhile, the algorithm learns the conditional network structure among outputs in an undirected Gaussian graphical model.

Different from existing techniques for multivariate regression/multi-task learning, **missoNet** allows missing values 
in the output data, gives reliable estimates without requiring users to have any prior knowledge for pre-processing 
them (e.g., imputation).


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
## generate data with overall missing probability 0.1, missing mechanism "MCAR"
sim.dat <- generateData(n = 300, p = 50, q = 20, rho = 0.1, missing.type = "MCAR")
tr <- 1:240
tst <- 241:300

## train the model using the first 240 samples
cv.obj <- cv.missoNet(X = sim.dat$X[tr, ], Y = sim.dat$Z[tr, ], kfold = 5)

## or train the model in parallel
library(snowfall)
cv.obj <- cv.missoNet(X = sim.dat$X[tr, ], Y = sim.dat$Z[tr, ], kfold = 5,
                      parallel = TRUE, cpus = min(parallel::detectCores()-1, 5)) 

## plot the cross-validation error heatmap
plot(cv.obj)

## estimated covariate coefficients that give the minimum CV error.
Bhat <- cv.obj$est.min$Beta

## make prediction
newy <- predict(cv.obj, newx = sim.dat$X[tst, ], s = "lambda.min")
```


## Learn more


## References

- [Cocolasso for high dimensional error-in-variables regression](https://arxiv.org/pdf/1510.07123.pdf)

- [High-dimensional regression with noisy and missing data: Provable guarantees with nonconvexity](https://arxiv.org/pdf/1109.3714.pdf)

