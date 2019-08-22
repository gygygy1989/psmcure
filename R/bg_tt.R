#' Generate Weibull Parameters from Life Table
#'
#' Generate Weibull Parameters from Life Table
#'
#' @param weipar weibull parameters
#' @param try the data set of survival time and survival probabilities
#' @return None
#'
#'
#'
lifett<-function(weipar,try){
  lambda<-exp(weipar[1])
  alpha<-exp(weipar[2])
  ll<-sum(abs(1-pweibull(try$time, lambda, alpha)-try$surv_mean_usa))
  return(ll)
}
