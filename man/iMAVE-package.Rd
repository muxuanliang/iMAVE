\name{iMAVE-package}
\alias{iMAVE-package}
\alias{iMAVE}
\docType{package}
\title{
\packageTitle{iMAVE}
}
\description{
\packageDescription{iMAVE}
}
\details{
This pakcage implements the iMAVE methods in Liang&Yu's paper (see reference). The iMAVE method achieves the dimension reduction for the model effect modification (intercation between a Treatment and many modifiers). It includes three functions, iMAVE, cv.iMAVEdim, and iMAVEInferAugment. The iMAVE implements iMAVE and iMAVE2 with a pre-specified dimension. The cv.iMAVEdim selects the optimal dimension by cross-validation. The iMAVEInferAugment provides the se estimation of the iMAVE or iMAVE2 when dimension is 1.
}
\author{
\packageAuthor{iMAVE}

Maintainer: \packageMaintainer{iMAVE}
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
#### sampling
nobs <- 500
nvars <- 10
gamma <- 0.1
tau <- 7
sigma <- 0.6

x <- matrix(rnorm(nobs * nvars), nobs, nvars)
beta1 <- c(1,1,1,1,1,1,1,1,1,1)
beta2 <- c(1,-1,1,-1,1,-1,1,-1,1,-1)

#### generate interaction and main effect ########
inter_eff <- tau * (pnorm(x %*% beta1) - 0.5) + tau * (pnorm(x %*% beta2) - 0.5)
main_eff <- (gamma * x %*% beta1)^2

#### generate noise and data #####################
epsilon <- rnorm(nobs, 0, sigma)

Tr <- rnorm(nobs)>0

y <- (Tr-0.5) * inter_eff + main_eff + epsilon

#### run iMAVEInferAugment ###########
fit <- iMAVEInferAugment(fit=NULL, x, y, Tr, pi = rep(0.5, nobs), d=1)
}
