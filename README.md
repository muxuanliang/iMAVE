# R Package: iMAVE

This pakcage implements the iMAVE methods in Liang&Yu's paper (see [reference](https://arxiv.org/abs/1804.05373)). The iMAVE method achieves the dimension reduction for the model effect modification (intercation between a Treatment and many modifiers). It includes three functions, iMAVE, cv.iMAVEdim, and iMAVEInferAugment. The iMAVE implements iMAVE and iMAVE2 with a pre-specified dimension. The cv.iMAVEdim selects the optimal dimension by cross-validation. The iMAVEInferAugment provides the se estimation of the iMAVE or iMAVE2 when dimension is 1.

# Details

Version: 1.0

Date: 2020-06-06

Author: Muxuan Liang <mliang@fredhutch.org>

Maintainer: Muxuan Liang <mliang@fredhutch.org>

Description: Dimensional reduction for the model effect modification. It utilizes MAVE type methods to reduce the dimension for the treatment-covariates interaction. When assumed no unmeasured confounding, it is connected to the individialized treatment effect.

License: GPL-3

Imports: Rcpp (>= 1.0.3),
         RcppEigen,
         mgcv,
         pracma

LinkingTo: Rcpp, RcppEigen

# Important Notice for Windows User

This R package involves C++ code. Before you compile from the source, please install the [Rtools40](https://cran.r-project.org/bin/windows/Rtools/) tool. 

# Example
```
#### Sampling
nobs <- 500

nvars <- 10

gamma <- 0.1

tau <- 7

sigma <- 0.6

x <- matrix(rnorm(nobs * nvars), nobs, nvars)

beta1 <- c(1,1,1,1,1,1,1,1,1,1)

beta2 <- c(1,-1,1,-1,1,-1,1,-1,1,-1)

#### Generate interaction and main effect ########
inter_eff <- tau * (pnorm(x %*% beta1) - 0.5) + tau * (pnorm(x %*% beta2) - 0.5)

main_eff <- (gamma * x %*% beta1)^2

#### Generate noise and data #####################
epsilon <- rnorm(nobs, 0, sigma)

Tr <- rnorm(nobs)>0

y <- (Tr-0.5) * inter_eff + main_eff + epsilon

#### Run iMAVE with dimension selection ###########
fit_select_dim <- cv.iMAVEdim(x, y, Tr, pi = rep(0.5, nobs))
fit <- iMAVE(x, y, Tr, pi = rep(0.5, nobs), d=fit_select_dim$opt_dim)
```
