#' Likelihood function for latency survival function in variance estimation
#'
#' Likelihood function for latency survival function in variance estimation
#'
#' @param beta The estimated parameters in latency survival function.
#' @param s The estimated baseline survival rate.
#' @param th  The estiamted total baseline hazard rate.
#' @param Status Right censoring status: 1 is death, 0 is censored.
#' @param X The covaraites for latency survival function.
#' @param w unobserved cure status indicator.
#' @param mort.h background hazard rate from life table.
#' @param mort.s background survival rate from life table.
#' @param ii the number of observation.
#'
#' @return None
#'
#'
#'
fish_func2<- function(beta,s,th,Status,X,w,mort.h,mort.s,ii){
  if (is.matrix(X)){
    temp2<- exp(beta%*%(X[ii,]))
  }else{
    temp2<- drop(exp(beta%*%X[ii]))
  }
  ss<- s[ii]^temp2
  thh<- th[ii]*temp2
  ll2<- Status[ii]*w[ii]*log(thh+(1e-10))+w[ii]*log(mort.s[ii]*ss+(1e-10))+Status[ii]*(1-w[ii])*log(mort.h[ii]+(1e-10))+(1-w[ii])*log(mort.s[ii]+(1e-10))
  return(ll2)
}
