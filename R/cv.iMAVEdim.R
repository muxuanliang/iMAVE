#' cv-MAVEdim for dimension selection in iMAVE or IMAVE2
#' 
#' @param x input matrix of dimension nobs \times nvars. Each raw is a observation, each column is a covariate
#' @param y numeric response
#' @param tr is a vector of binary value representing two treatment, 0 or 1. 
#' @param pi is the propensity score p(tr=1|X)
#' @param nfolds number of the folds in the cross-validation
#' @param parallel implies whether use parallel computing in the cross-validation
#' @return An object with S3 class "iMAVE"
#' @import Rcpp
#' 
#' @export
#' 
cv.iMAVEdim <- function(x, y, tr,
                      dims = c(1:ncol(x)),
                      pi = 0.5,
                      nfolds = 5,
                      parallel = F,
                      ...)
{
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
  
  
  
  numdim <- length(dims)
  foldid=sample(rep(seq(nfolds),length=nobs))
  if(nfolds<3)stop("nfolds must be bigger than 3; nfolds=5 recommended")
  outlist <- list()
  for (numdim_select in 1:numdim){
    outlist[[numdim_select]] <- list()
  }
  
  if (parallel && require(foreach)) {
    library(doParallel)
    n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
  }
  
  for (numdim_select in 1:numdim){
    if (parallel && require(foreach)) {
      outlist[[numdim_select]] = foreach (i=seq(nfolds), .packages=c("iMAVE")) %dopar% {
        which=foldid==i
        iMAVE(x[!which,,drop=FALSE], y[!which], tr[!which], pi = pi[!which], d = dims[numdim_select], ...)
      }
    }else{
      for(i in seq(nfolds)){
        which=foldid==i
        outlist[[numdim_select]][[i]]=iMAVE(x[!which,,drop=FALSE], y[!which], tr[!which], pi = pi[!which], d = dims[numdim_select], ...)
      }
    }
  }
  
  if (parallel && require(foreach)) {
    stopCluster(cl)
  }
  
  predmat=matrix(NA,length(y),numdim)
  cvraw <- array(0,c(nfolds,numdim))
  for(j in seq(numdim)){
    fitobj=outlist[[j]]
    cvraw[,j] <- loss.score(nfolds, foldid, fitobj, x, y, tr, pi)
  }
  
  cvm=apply(cvraw,2,mean,na.rm=TRUE)
  cvsd=apply(cvraw,2,sd,na.rm=TRUE)
  opt.index <- max(seq(1, numdim, by = 1)[cvm==cvm[which.min(cvm)]])
  
  iMAVE.call=match.call(expand.dots=TRUE)
  which=match(c("nfolds"),names(iMAVE.call),F)
  
  if(any(which))iMAVE.call = iMAVE.call[-which]
  
  iMAVE.call[[1]] = as.name("iMAVE")
  iMAVE.object = iMAVE(x, y, tr, pi = pi, d = dims[opt.index], ...)
  iMAVE.object$call=iMAVE.call
  
  out=list(cvm=cvm,cvsd=cvsd,cvup=cvm+cvsd,
           cvlo=cvm-cvsd, cvraw = cvraw, iMAVE.fit=iMAVE.object, beta = iMAVE.object$beta)
  class(out)="cv.iMAVE"
  out
}

loss.score <- function(nfolds, foldid, fitobj, x, y, tr, pi){
  
  loss <- array(0,c(nfolds, 1))
  for (i in 1:nfolds){
    which = foldid==i
    x.train <- x[!which,,drop=F]
    y.train <- y[!which]
    tr.train <- tr[!which]
    pi.train <- pi[!which]
    x.test <- x[which,,drop=F]
    y.test <- y[which]
    tr.test <- tr[which]
    pi.test <- pi[which]
    
    beta <- fitobj[[i]]$beta
    
    if(anyNA(beta)){
      loss[i] <- NA
    } else {
      pred.link <- x.test %*% beta
      train.link <- x.train %*% beta
      
      y.train.weighted <- (tr.train-0.5) * y.train / (pi.train * tr.train + (1-pi.train) * (1-tr.train))
      
      pred <- ks(train.link,y.train.weighted, pred.link)
      yy.test <- (tr.test - 0.5) * y.test/(pi.test * tr.test + (1-pi.test) * (1-tr.test))
      loss[i] <- mean((yy.test-pred)^2)
    }
  }
  loss
}

ks <- function(x, y, x.test){
  nobs <- nrow(x)
  nvars <- ncol(x)
  hopt <- (4/(nvars+2))^(1/(nvars+4)) * (nobs^(-1/(nvars+4)))
  wm <- function(t){
    if (ncol(x)==1){
      weight <- exp(-0.5 * (t-x)^2/(hopt^2)) * hopt
    } else {
      weight <- apply(x,1,function(s){exp(-0.5 * sum((t-s)^2)/(hopt^2)) * hopt^(ncol(x))})
    }
    weighted.mean(y, weight)
  }
  y.test<- array(0,c(nrow(x.test),1))
  for (index in 1:nrow(x.test)){
    if (ncol(x.test)==1){
      y.test[index] <- wm(x.test[index])
    } else {
      y.test[index] <- wm(x.test[index,])
    }
  }
  y.test
}