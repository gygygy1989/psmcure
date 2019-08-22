#' Complete Likelihood Function in Parametric Method
#'
#' Complete Likelihood Function in Parametric Method
#'
#' @param par all the parameters in the model.
#' @param Time Survival time from data set.
#' @param Status Right censoring status: 1 is death, 0 is censored.
#' @param Z The covariates included in the cure linke function.
#' @param X The covaraites for latency survival function.
#' @param dis If "parametric" is specified for method, then dis must be provided, such as "exp","weibull","llogis","gamma","gompertz","lnorm" and "gengamma".
#' @param w unobserved cure status indicator.
#' @param link The type of cure rate link function, including "logit","probit", and "cloglog".
#' @param mort.h background hazard rate.
#' @param mort.s background survival rate.
#' @param ii the number of observation.
#'
#' @return None
#'
#'
#'
com_like<- function(par,Time,Status,Z,X,dis,w,link,mort.h,mort.s,ii){
  b<- par[1:ncol(Z)]
  if(is.matrix(X)){
    beta<- par[(ncol(Z)+1):(ncol(Z)+ncol(X))]
    dpar<- par[-c(1:(ncol(Z)+ncol(X)))]
  }else{beta<- par[(ncol(Z)+1):(ncol(Z)+NCOL(X))]
  dpar<- par[-c(1:(ncol(Z)+NCOL(X)))]}


  temp1<- exp(b%*%Z[ii,])
  if(link=="logit"){rate<- temp1/(1+temp1)}
  if(link=="probit") {rate<- pnorm(log(temp1))}
  if(link=="cloglog"){rate<- 1-exp(-temp1)}

  if (is.matrix(X)){
    temp2<- exp(beta%*%X[ii,])
  }else{temp2<- exp(beta%*%X[ii])}

  base_surv<-psurv(Time=Time[ii],dis=dis,dpar=dpar)$basesurv
  base_haz<- psurv(Time=Time[ii],dis=dis,dpar=dpar)$basehaz
  ss<- base_surv^(temp2)
  hh<- base_haz*(temp2)
  part1<- w[ii]*Status[ii]*log(mort.h[ii]+hh)+w[ii]*log(ss)
  part2<- mort.h[ii]^Status[ii]*mort.s[ii]
  lpart2<- ifelse(part2==0,0,log(part2))
  logl<- w[ii]*log(rate)+(1-w[ii])*log(1-rate)+(1-w[ii])*lpart2+ifelse(ss==0,0,part1)+w[ii]*log(mort.s[ii])
  return(logl)
}
