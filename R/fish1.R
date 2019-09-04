#' Likelihood function for cure link function in variance estimation.
#'
#' Likelihood function for cure link function in variance estimation.
#'
#' @param b parameters in cure link function.
#' @param Z The covariates included in the cure link function.
#' @param w unobserved cure status indicator.
#' @param link The type of cure rate link function, including "logit","probit", and "cloglog".
#' @param ii the number of observation.
#' @importFrom stats pnorm
#' @return None
#'
#'
#'
fish_func1<- function(b,Z,w,link,ii){
  temp1<- exp(b%*%(Z[ii,]))
  if(link=="logit") {rate<- as.vector(temp1/(1+temp1)) }
  if(link=="probit"){rate<- as.vector(pnorm(log(temp1)))}
  if(link=="cloglog"){rate<- as.vector(1-exp(-temp1))}
  ll1<- w[ii]*log(rate)+(1-w[ii])*log(1-rate)
}
