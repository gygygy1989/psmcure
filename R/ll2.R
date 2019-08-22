#' Likelihood function for latency survival function in EM algorithm
#'
#' Likelihood function for latency survival function in EM algorithm
#'
#' @param trybeta The estimated parameters in latency survival function.
#' @param Time Survival time from data set.
#' @param Status Right censoring status: 1 is death, 0 is censored.
#' @param age the variable of age
#' @param X The covaraites for latency survival function.
#' @param w unobserved cure status indicator.
#' @param mort.s background survival rate from life table.
#' @param mort.h background hazard rate from life table.
#'
#' @return None
#'
#'
#'
ll2<-function(trybeta,Time,Status,age,X,w,mort.s,mort.h){
  trybase<- surv_func(Time,Status,age,X,trybeta,w,mort.s,mort.h)
  trys<-  trybase$basesurv
  tryth<-  trybase$haz
  exb<- drop(exp(trybeta%*%t(X)))
  tryss<- drop(trys^exb)
  trythh<- drop(tryth*exb)
  -sum(Status*w*log(trythh+1e-10)+w*log(mort.s*tryss+1e-10))
}
