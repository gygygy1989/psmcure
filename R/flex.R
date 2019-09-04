#' Flexible Parametric EM Algorithm Main Function.
#'
#' Flexible Parametric EM Algorithm Main Function.
#'
#' @param Time Survival time from data set.
#' @param Status Right censoring status: 1 is death, 0 is censored.
#' @param Z The covariates for cure link function.
#' @param X The covaraites for latency survival function.
#' @param Mbigs M spline basis.
#' @param Ibigs I spline basis.
#' @param r Initial value of r.
#' @param b Initial value of b.
#' @param beta Initial value of beta.
#' @param w Initial value of w.
#' @param mort.s background survival rate from life table.
#' @param mort.h background hazard rate from life table.
#' @param link The type of cure rate link function, including "logit","probit", and "cloglog".
#' @param emmax The maximum iteration is 100.
#' @param eps The stopping criteria for convergence is 1e-7.
#' @importFrom stats glm coef pnorm optim
#' @importFrom numDeriv grad
#' @importFrom MASS ginv
#' @return None
#'
#'
#' @export
#'
flex.spline<- function(Time,Status,Z,X,Mbigs,Ibigs,r,b,beta,w,mort.s,mort.h,link="logit",emmax=100,eps=1e-10){

  convergence=1000;i=1
  while(convergence>eps & i<emmax){

    if(link=="logit"){uncure=as.vector(exp(b%*%t(Z))/(1+exp(b%*%t(Z))))}
    if(link=="probit"){uncure=as.vector(pnorm(b%*%t(Z)))}
    if(link=="cloglog"){uncure=as.vector(1-exp(-exp(b%*%t(Z))))}

    h<- as.numeric(r%*%Mbigs)
    cbase<- as.numeric(r%*%Ibigs)
    s<- exp(-cbase)
    exb<- drop(exp(beta%*%t(X)))

    su<- s^exb
    hu<- h*exb

    ##E step
    w=Status*(uncure*(su)*(mort.h+hu))/((1-uncure)*mort.h+uncure*(su)*(mort.h+hu))+
      (1-Status)*(uncure*su)/((1-uncure)+(uncure*su))
    ###########

    ##M step
    update_b = optim(par=b,fn=ll1r,Z=Z,w=w,link=link)$par
    update_par = optim(par=c(beta,log(r+1e-7)),fn=ll2r,Time=Time,Status=Status,X=X,Mbigs=Mbigs,Ibigs=Ibigs,w=w,mort.s=mort.s,mort.h=mort.h)$par
    update_beta=update_par[1:length(beta)]
    update_r=exp(update_par[-c(1:length(beta))])
    #####################################

    convergence=sum(c((update_b-b)^2),(c(update_beta-beta)^2))

    if(convergence==Inf || is.na(convergence)) {
      break
    }else {
      b<-update_b
      beta<- update_beta
      r<- update_r
    }
    i<- i+1
    #cat("i=",i,"cov=",convergence)
  }

  ##### variance
  est0= c(b,beta)
  est0r=c(b,beta,r)
  thessian<-diag(rep(0,length(est0r)))
  n<-length(Time)
  for(ii in 1: n){
    tempg<-grad(like.com.one,x=est0r,link=link,X=X,Z=Z,Mbigs=Mbigs,Ibigs=Ibigs,Status=Status,w=w,mort.s=mort.s,mort.h=mort.h,ii=ii)
    tempm<-tempg%*%t(tempg)
    thessian<-thessian+tempm
  }
  se0r<- sqrt(diag(ginv(thessian)))
  se0<- se0r[1:(length(b)+length(beta))]

  h<- as.numeric(r%*%Mbigs) #h0
  cbase<- as.numeric(r%*%Ibigs) #H0
  s<- exp(-cbase)  #s0

  flex.out<- list(b=b,latencyfit=beta,r=r,se=se0,Uncureprob=uncure,w=w,conv=convergence,s=s,h=h)
}

