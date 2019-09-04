#' Second Likelihood function for latency survival function in Flexible EM algorithm.
#'
#' Second Likelihood function for latency survival function in Flexible EM algorithm.
#'
#' @param x The parameters to be updated.
#' @param Time Survival time from data set.
#' @param Status Right censoring status: 1 is death, 0 is censored.
#' @param X The covaraites for latency survival function.
#' @param Mbigs M spline basis.
#' @param Ibigs I spline basis.
#' @param w unobserved cure status indicator.
#' @param mort.s background survival rate from life table.
#' @param mort.h background hazard rate from life table.
#' @return None
#'
#'
#'
#'
ll2r<- function(x,Time,Status,X,Mbigs,Ibigs,w,mort.s,mort.h){

  beta<- x[1:ncol(X)]
  exb<- drop(exp(beta%*%t(X)))

  r<- exp(x[-c(1:ncol(X))])
  h<- as.numeric(r%*%Mbigs) #h0
  cbase<- as.numeric(r%*%Ibigs) #H0
  s<- exp(-cbase)  #s0

  -sum(w*Status*log(mort.h+h*exb+(1e-7))-w*exb*cbase)
}
