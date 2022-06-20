---
title: "README"
output: html_document
---

# missoNet: a package for multivariate regression and conditional network estimation with missing outputs.

The package **missoNet** enables users to estimate covariate effects by leveraging the information of multiple correlated outputs. 
Meanwhile, the algorithm learns the network structure of outputs in an undirected Gaussian graphical model.

Different from existing multivariate methods, **missoNet** allows the output data to contain missing values, and does not require users to preprocess them (e.g., imputation).


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
cv.obj <- cv.missoNet(X = sim.dat$X, Y = sim.dat$Z, kfold = 5) ## run the analysis

plot(cv.obj)  ## plot the cross-validation error surface

newy <- predict(cv.obj, newx = sim.dat$X, s = "lambda.min")  ## make prediction
```


## Learn more


## References
