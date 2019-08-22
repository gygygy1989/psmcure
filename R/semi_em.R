#' Semiparametric EM algorithm to estimate parameters
#'
#' Semiparametric EM algorithm to estimate parameters
#'
#' @param Time Survival time from data set.
#' @param Status Right censoring status: 1 is death, 0 is censored.
#' @param age The varaible of age.
#' @param Z The covariates for cure link function.
#' @param X The covaraites for latency survival function.
#' @param mort.s background survival rate from life table.
#' @param mort.h background hazard rate from life table.
#' @param link The type of cure rate link function, including "logit","probit", and "cloglog".
#' @param emmax The maximum iteration is 100.
#' @param eps The stopping criteria for convergence is 1e-7.
#'
#' @return None
#'
#'
#' @export
emback<-function(Time,Status,age,Z,X,mort.s,mort.h,link="logit",emmax=100,eps=1e-7){
  Time=as.numeric(Time==0)*(1e-10)+as.numeric(Time!=0)*Time

  ## initial values
  w= Status
  n=length(Status)
  b=glm(w~Z[,-1],family=quasibinomial(link=link))$coef
  beta=coxph(Surv(Time,Status)~X,method="breslow")$coef
  ############################
  ############################

  convergence=1000;i=1

  while(convergence>eps & i<emmax){

    if(link=="logit"){uncure=as.vector(exp(b%*%t(Z))/(1+exp(b%*%t(Z))))}
    if(link=="probit"){uncure=as.vector(pnorm(b%*%t(Z)))}
    if(link=="cloglog"){uncure=as.vector(1-exp(-exp(b%*%t(Z))))}

    s<- surv_func(Time,Status,age,X,beta,w,mort.s,mort.h)$basesurv
    th<- surv_func(Time,Status,age,X,beta,w,mort.s,mort.h)$haz


    su<- drop(s^exp(beta%*%t(X)))
    thazard<- drop(th*exp(beta%*%t(X)))

    ##E step
    w=Status*(uncure*su*thazard/((1-uncure)*mort.h+uncure*su*thazard))+(1-Status)*((uncure*su)/((1-uncure)+(uncure)*su))
    ###########

    ##M step
    #update b
    update_b=optim(par=b,fn=ll1,Z=Z,w=w,link=link)$par

    #update beta
    if(length(beta)==1){
      update_beta = optimise(f=ll2,interval=c((beta-5),(beta+5)),Time=Time,Status=Status,age=age,X=X,w=w,mort.s=mort.s,mort.h=mort.h)$minimum
    }else if (length(beta)>1){
      update_beta = optim(par=beta,fn=ll2,Time=Time,Status=Status,age=age,X=X,w=w,mort.s=mort.s,mort.h=mort.h)$par
    }

    convergence=sum(c((update_b-b)^2),c((update_beta-beta)^2))

    if(convergence==Inf || is.na(convergence)) {
      break
    }else {
      b<-update_b
      beta<- update_beta
    }
    i<- i+1
  }
  ######################################

  #variance
  se<- per_fun(b=b,Z=Z,link=link,Status=Status,w=w,beta=beta,s=s,th=th,X=X,mort.h=mort.h,mort.s=mort.s)
  bse<- se[1:length(b)]
  betase<- se[(length(b)+1):(length(b)+length(beta))]

  emback <- list(b=b,latencyfit=beta,Uncureprob=uncure,s=s,th=th,bse=bse,betase=betase,w=w)
}


