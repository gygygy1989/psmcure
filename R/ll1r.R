#' First Likelihood Function to Update b in Flexible Parametric.
#'
#' First Likelihood Function to Update b in Flexible Parametric.
#'
#' @param tryb The value of b to be updated.
#' @param Z The covariates for cure link function.
#' @param w Unobserved cure status indicator.
#' @param link The type of cure rate link function, including "logit","probit", and "cloglog".
#' @importFrom stats pnorm
#' @return None
#'
#'
#'
#'
ll1r<-function(tryb,Z,w,link){
  if(link=="logit"){uncure<- as.vector(exp(tryb%*%t(Z))/(1+exp(tryb%*%t(Z))))}
  if(link=="probit"){uncure<- as.vector(pnorm(tryb%*%t(Z)))}
  if(link=="cloglog"){uncure<- as.vector(1-exp(-exp(tryb%*%t(Z))))}
  -sum(w*log(uncure+1e-10)+(1-w)*log(1-uncure+1e-10))
}
