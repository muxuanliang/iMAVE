#' iMAVEInferAugment for the infetence of the estimated index
#' 
#' @param fit can be a iMAVE object to be used for inference. Defaults run iMAVE internally.
#' @param x input matrix of dimension nobs x nvars. Each raw is a observation, each column is a covariate
#' @param y numeric response
#' @param tr is a vector of binary value representing two treatment, 0 or 1. 
#' @param pi is the propensity score (success)
#' @param efficient imdicates whther use efficient score when the condition variance of the outcome is a constant
#' @param ... can be other parameters in iMAVE function
#' @return A list
#' \describe{
#' \item{coef_one_step}{Estimated direction after one-step estimation.}
#' \item{se}{Estimated se of \code{coef_one_step}.}
#' }
#' @import Rcpp
#'
#' @export
#' 
iMAVEInferAugment <- function(fit = NULL, x, y, tr, pi = 0.5, efficient = FALSE,  ...){
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

