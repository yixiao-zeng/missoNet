---
title: "README"
output: html_document
---

# missoNet: a R package for multivariate regression and conditional network estimation with missing outputs.

The package **missoNet** enables users to estimate covariate effects by leveraging the information of multiple correlated outputs. 
Meanwhile, the algorithm learns the network structure among outputs in an undirected Gaussian graphical model.

Different from existing multivariate regression/multi-task learning techniques, **missoNet** allows the output data to contain missing values, and does not require users to have any prior knowledge for preprocessing them (e.g., imputation).


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
sim.dat <- generateData(n = 300, p = 50, q = 20, rho = 0.1, missing.type = "MCAR")  ## generate data
cv.obj <- cv.missoNet(X = sim.dat$X, Y = sim.dat$Z, kfold = 5)  ## train the model

plot(cv.obj)  ## plot the cross-validation error surface

Bhat <- cv.obj$est.min$Beta  ## estimated covariate coefficients that give the minimum CV error.

newy <- predict(cv.obj, newx = sim.dat$X, s = "lambda.min")  ## make prediction
```


## Learn more


## References

- [Cocolasso for high dimensional error-in-variables regression](https://arxiv.org/pdf/1510.07123.pdf)

- [High-dimensional regression with noisy and missing data: Provable guarantees with nonconvexity](https://arxiv.org/pdf/1109.3714.pdf)

