#' Likelihood function for cure link function in EM algorithm
#'
#' Likelihood function for cure link function in EM algorithm
#'
#' @param tryb parameters in cure link function.
#' @param Z The covariates included in the cure link function.
#' @param w unobserved cure status indicator.
#' @param link The type of cure rate link function, including "logit","probit", and "cloglog".
#'
#' @return None
#'
#'
#'
ll1<-function(tryb,Z,w,link){
  if(link=="logit"){uncure<- as.vector(exp(tryb%*%t(Z))/(1+exp(tryb%*%t(Z))))}
  if(link=="probit"){uncure<- as.vector(pnorm(tryb%*%t(Z)))}
  if(link=="cloglog"){uncure<- as.vector(1-exp(-exp(tryb%*%t(Z))))}
  -sum(w*log(uncure)+(1-w)*log(1-uncure))
}
