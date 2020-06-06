#' iMAVEInferAugment for the infetence of the estimated index
#' 
#' @param x input matrix of dimension nobs \times nvars. Each raw is a observation, each column is a covariate
#' @param y numeric response
#' @param tr is a vector of binary value representing two treatment, 0 or 1. 
#' @param pi is the propensity score (success)
#' @param efficient imdicates whther use efficient score when the condition variance of the outcome is a constant
#' @import Rcpp
#'
#' @export
#' 
iMAVEInferAugment <- function(fit = NULL, x, y, tr, pi = 0.5, efficient = FALSE, ...){
   n <- length(y)
   sampleSplitIndex <- (rnorm(n)>0)
   
   if (is.null(fit)){
      fit <- iMAVE(x, y, tr, pi, d=1, ...)
   }
   
   mainEffect <- NULL
   if(class(x)=="iMAVE2"){
      mainEffect <- fit$mainEffect
   }
   fit1 <- iMAVEInfer(fit, x, y, tr, pi, mainEffect = mainEffect, sampleSplitIndex = sampleSplitIndex, efficient = efficient)
   fit2 <- iMAVEInfer(fit, x, y, tr, pi, mainEffect = mainEffect, sampleSplitIndex = !sampleSplitIndex, efficient = efficient)
   res <- list(coef=(fit1$coef + fit2$coef)/2, coef_one_step=(fit1$coef_one_step+fit2$coef_one_step)/2, se=sqrt((fit1$se^2 + fit2$se^2)/2)/sqrt(n))
   res
}

iMAVEInfer <- function(fit, x, y, tr, pi, mainEffect = NULL, sampleSplitIndex = NULL, efficient = FALSE)
{
 d <- dim(fit$beta)[2]
 p <- dim(fit$beta)[1]
 n <- length(y)
 if (dim(fit$beta)[2] > 1){
   coef <- fit$beta %*% solve(fit$beta[1:d,1:d])
 } else {
   coef <- fit$beta/fit$beta[1]
 }
 
 ks <- function(xx, yy, xx.test){
   nobs <- nrow(xx)
   nvars <- ncol(xx)
   hopt <- (4/(nvars+2))^(1/(nvars+4)) * (nobs^(-1/(nvars+4)))
   wm <- function(t){
     if (ncol(xx)==1){
       weight <- exp(-0.5 * (as.numeric(t)-xx)^2/(hopt^2)) * hopt
     } else {
       weight <- apply(xx,1,function(s){exp(-0.5 * sum((t-s)^2)/(hopt^2)) * hopt^(ncol(xx))})
     }
     weighted.mean(yy, weight)
   }
   if (nrow(as.matrix(xx.test))==1) {
     yy.test <- wm(xx.test)
   } else {
     yy.test<- array(0,c(nrow(xx.test),1))
     for (index in 1:nrow(xx.test)){
       if (ncol(xx.test)==1){
         yy.test[index] <- wm(xx.test[index])
       } else {
         yy.test[index] <- wm(xx.test[index,])
       }
     }
   }
   yy.test
 }
 
 iMAVEpredict <- function(x.link){
   x.link <- rbind(x.link,x.link)
   predict(shadow_model, newdata = list(predictor = x.link),type='response')[1]
 }
 
 exp_gam <- function(x, y){
    shadow_data <- list(predictor = x, response = y)
    shadow_model <- lm(response~predictor, data = shadow_data)
    as.matrix(predict(shadow_model, newdata = list(predictor = link),type='response'))
 }
 
 xNormalized <- function(x.test){
   pred.link <- x.test %*% coef
   train.link <- x %*% coef
   x.test.normalized <- x.test * 0
   for (i in 1:length(x.test)){
     x.test.normalized[i] <- ks(train.link,x[,i], pred.link)
   }
   x.test #- x.test.normalized
 }
 
 iMAVEgradient <- function(x.link){
   pracma::grad(iMAVEpredict, x.link)
 }
 
 # set sampleSplitIndex
 if(is.null(sampleSplitIndex)){
    sampleSplitIndex <- (rnorm(n) > 0)
 }
 
 link <- x %*% coef
 V <- rep(1, times=n)
 EV <- rep(1, times=n)
 if (efficient){
    V <- pi * (1-pi)
    EV <- exp_gam(link, V)
 }
 
 shadow_data <- list(predictor = link[sampleSplitIndex,,drop=FALSE], response = 2 * (tr[sampleSplitIndex]-0.5) * y[sampleSplitIndex]/(pi[sampleSplitIndex]*tr[sampleSplitIndex]+(1-pi[sampleSplitIndex])*(1-tr[sampleSplitIndex])))
 if(nrow(fit$ab)>=n){
   shadow_data$response <- fit$ab[sampleSplitIndex,1]
 }
 shadow_model <- mgcv::gam(response~predictor, data = shadow_data)
 g.value <- as.matrix(predict(shadow_model, newdata = list(predictor = link),type='response'))
 g.gradient <- apply(link, 1, iMAVEgradient)
 EX <- apply(x, 2, function(t){exp_gam(link, V*t)})
 
 if (is.null(mainEffect)){
   residual <- y-g.value*(tr-0.5)
 } else {
   shadow_data <- list(predictor = link[sampleSplitIndex,], response = 0.5*y[sampleSplitIndex]/(pi[sampleSplitIndex]*tr[sampleSplitIndex]+(1-pi[sampleSplitIndex])*(1-tr[sampleSplitIndex])))
   shadow_model <- mgcv::gam(response~predictor, data = shadow_data)
   mainEffect <- predict(shadow_model, newdata = list(predictor = link),type='response')
   residual<-(y-g.value*(tr-0.5)-as.matrix(mainEffect))
 }
 A_tmp <- pracma::zeros(d*(p-d))
 S_tmp <- pracma::zeros(d*(p-d), 1)
 B_tmp <- pracma::zeros(d*(p-d))
 for (iter in (c(1:n)[!sampleSplitIndex])){
   if (d==1){
     A.left_tmp <- t(kronecker(t(g.gradient[iter]), x[iter,-(1:d)]-EX[iter,-(1:d)]/EV[iter]))
   } else {
     A.left_tmp <- matrix(t(kronecker(t(g.gradient[,iter]), x[iter,-(1:d)]-EX[iter,-(1:d)]/EV[iter])), c(1,d*(p-d)))
   }
   A_tmp <- 1/4 * t(A.left_tmp) %*% (A.left_tmp)/(pi[iter]*tr[iter]+(1-pi[iter])*(1-tr[iter]))* (V[iter]) + A_tmp
   S_tmp <- t(A.left_tmp) * (tr[iter]-0.5) * residual[iter] /(pi[iter]*tr[iter]+(1-pi[iter])*(1-tr[iter])) * V[iter] + S_tmp
   B_tmp <- 1/4 * t(A.left_tmp) %*% (A.left_tmp)/(pi[iter]*tr[iter]+(1-pi[iter])*(1-tr[iter]))^2 * (residual[iter])^2 * (V[iter])^2 + B_tmp
 } 
 A <- A_tmp/sum(!sampleSplitIndex)
 S <- S_tmp/sum(!sampleSplitIndex)
 B <- B_tmp/sum(!sampleSplitIndex)
 covMatrix <- pracma::inv(A) %*% B %*% pracma::inv(A)
 coef_one_step <- rbind(pracma::eye(d), coef[-d,] + pracma::inv(A) %*% S)
 covMatrix_se <- rbind(pracma::zeros(d),matrix(sqrt(diag(covMatrix)),c(p-d,p-d)))
 res <- list(coef=coef, coef_one_step=coef_one_step, covMatrix=covMatrix,se=covMatrix_se)
 res
}
