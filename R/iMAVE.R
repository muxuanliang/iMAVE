#' iMAVE for model effect modification
#' 
#' @param x input matrix of dimension nobs x nvars. Each raw is a observation, each column is a covariate
#' @param y numeric response
#' @param tr is a vector of binary value representing two treatment, 0 or 1. 
#' @param d number of the dimension to reduce
#' @param pi is the propensity score p(tr=1|X)
#' @param method  'zero' or 'all' decides the local linear expansion only near zero or all data point. 'all' is default for imave and imave2.
#' @param monoLink indicates wehter the link function is assumed to be monotone. False by default.
#' @param mainEffect mainEffect is TRUE if iMAVE2; FALSE if iMAVE
#' @param initial selects the methds to obtain initial values. Default option "direct" start from 0; "zero" solves an imave with method='zero'
#' @param tol is the convergence tolerence
#' @param maxit is the maxmum iteration
#' @param normalizeWeight indicates whether the weight for kernels is normalized. Default is TRUE.
#' @param constraint indicates the constaint on B; by default, Grassmann Manifold
#' @return An object with S3 class "iMAVE"
#' @import Rcpp
#' 
#' @export
#' 
iMAVE <- function(x, y, tr,
                   d = 1,
                   pi = 0.5,
                   method = "all",
                   monoLink = FALSE,
                   mainEffect = TRUE,
                   initial = "direct",
                   tol = 1e-5,
                   maxit = 100L,
                   normalizeWeight = TRUE,
                   constraint = "Grassmann")
{
  stopifnot(!((constraint == "norm") & (d>1)))
  method <- match.arg(method)
  
  if(missing(initial)|method == "zero")
  {
    initial <- "direct"
  }
  stopifnot(!(method == "zero" & initial == "zero"))
  initial <- match.arg(initial)
  
  if(monoLink == F)
  {
    const = 0
  } else {
    const = 1
  }
  
  nvars <- ncol(x)
  nobs <- nrow(x)
  y <- drop(y)
  dimy <- dim(y)
  leny <- ifelse(is.null(dimy), length(y), dimy[1])
  stopifnot(leny == nobs)
  
  if(length(pi)==1){
    pi = rep(pi, nobs)
  } else if(length(pi) < nobs){
    stop("Perpensity score need to be constant or a vector with length of observation")
  }
  
  stopifnot(all((as.numeric(tr) == 1) + (as.numeric(tr) == 0)))
  tr <- as.numeric(tr)
  
  opts <- list(tol = tol,
               maxit = maxit,
               constraint = const,
               normalizeWeight = normalizeWeight)
  
  fit <- .Call("iMAVE_cpp",
                 x_        = x,
                 y_        = y,
                 tr_       = tr,
                 d_        = d,
                 pi_       = pi,
                 init_     = initial,
                 method_   = method,
                 opts_     = opts, 
                 PACKAGE   = "iMAVE"
                 )
  if (d ==1){
    fit$beta <- switch(constraint, "Grassmann"=fit$beta/fit$beta[1], "norm"=fit$beta/sqrt(sum(fit$beta^2)))
  } else {
    fit$beta <- fit$beta %*% solve(fit$beta[1:d,1:d])
  }
  
  if (mainEffect){
      predict.link <- x %*% fit$beta
      predict.trEffect <- 2*ks(predict.link, y * (tr-0.5)/(pi*tr+(1-tr)*(1-pi)), predict.link)
      predict.epsilon <- y - (tr-0.5) * predict.trEffect
      predict.main <- ks(x, predict.epsilon/(pi^2*tr+(1-tr)*(1-pi)^2), x)
      predict.main.key <- 1/ks(x, 1/(pi^2*tr+(1-tr)*(1-pi)^2), x)*predict.main
      y.nomain <- y - predict.main.key
   
      fit <- .Call("iMAVE_cpp",
                   x_        = x,
                   y_        = y.nomain,
                   tr_       = tr,
                   d_        = d,
                   pi_       = pi,
                   init_     = initial,
                   method_   = method,
                   opts_     = opts, 
                   PACKAGE   = "iMAVE"
      )
      if (d ==1){
        fit$beta <- switch(constraint, "Grassmann"=fit$beta/fit$beta[1], "norm"=fit$beta/sqrt(sum(fit$beta^2)))
        fit$mainEffect <- y-y.nomain
      } else {
        fit$beta <- fit$beta %*% solve(fit$beta[1:d,1:d])
        fit$mainEffect <- y-y.nomain
      }
  }
  
  # if does not converge
  if (anyNA(fit$beta)){
    fit$beta <- array(0, c(nvars, d))
  }
  
  class2 <- switch (mainEffect,
    "TRUE" = "2", "FALSE" = ""
  )
  
  class(fit) <- paste0("iMAVE.", class2)
  fit
  }


  
  