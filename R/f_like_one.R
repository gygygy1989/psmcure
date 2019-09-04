#' Complete Likelihood Function For Each Subject in Flexible Parametric.
#'
#' Complete Likelihood Function For Each Subject in Flexible Parametric.
#'
#' @param est0r The parameters to be updated.
#' @param link Survival time from data set.
#' @param X The covariates for latency survival function.
#' @param Z The covariates fir cure function.
#' @param Mbigs M spline basis.
#' @param Ibigs I spline basis.
#' @param Status Right censoring status: 1 is death, 0 is censored.
#' @param w unobserved cure status indicator.
#' @param mort.s background survival rate from life table.
#' @param mort.h background hazard rate from life table.
#' @param ii The indicator for subject.
#' @return None
#'
#'
#'
#'
#'
like.com.one<- function(est0r,link,X,Z,Mbigs,Ibigs,Status,w,mort.s,mort.h,ii){
  b=est0r[1:ncol(Z)]
  beta=est0r[(ncol(Z)+1):(ncol(Z)+ncol(X))]
  r=est0r[-c(1:(ncol(Z)+ncol(X)))]

  if(link=="logit"){uncure=as.vector(exp(b%*%t(Z))/(1+exp(b%*%t(Z))))}
  if(link=="probit"){uncure=as.vector(pnorm(b%*%t(Z)))}
  if(link=="cloglog"){uncure=as.vector(1-exp(-exp(b%*%t(Z))))}

  h<- as.numeric(r%*%Mbigs)
  cbase<- as.numeric(r%*%Ibigs)
  s<- exp(-cbase)

  exb<- drop(exp(beta%*%t(X)))

  su<- s^exb
  hu<- h*exb

  logl<- w[ii]*log(uncure[ii])+(1-w[ii])*log(1-uncure[ii])+w[ii]*Status[ii]*log(mort.h[ii]+hu[ii])+w[ii]*log(mort.s[ii])+
    w[ii]*log(su[ii]+1e-7)+(1-w[ii])*Status[ii]*log(mort.h[ii])+(1-w[ii])*log(mort.s[ii])

  return(logl)
}
