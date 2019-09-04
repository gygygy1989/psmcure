#' Parametric EM Algorithm.
#'
#' Parametric EM Algorithm.
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
#' @param dis The parametric distribution such as "exp","weibull","llogis","gamma","gompertz","lnorm" and "gengamma".
#' @importFrom stats glm coef pnorm optim
#' @importFrom flexsurv flexsurvreg
#' @importFrom numDeriv grad
#' @importFrom MASS ginv
#' @return None
#'
#'
#' @export
emp<-function(Time,Status,age,Z,X,mort.s,mort.h,link="logit",emmax=100,eps=1e-7,dis=dis){

  if(is.null(dis)){stop("Please specify distribution!",call. = FALSE)
                   }
  Time=as.numeric(Time==0)*(1e-10)+as.numeric(Time!=0)*Time

  #####initial
  n=length(Status)
  b=glm((Status)~Z[,-1],family=quasibinomial(link=link))$coef

  ndpar=2
  if(dis=="exp") ndpar=1
  if(dis=="gengamma") ndpar=3

  betafit=flexsurvreg(Surv(Time,Status)~X,dist=dis)
  dpar=betafit$res[,1][1:ndpar]
  lgdpar<-betafit$res.t[,1][1:ndpar]   #log of part of parameters
  beta<-coef(betafit)[-c(1:ndpar)]


  #####start update
  i=1;convergence=1000
  while(convergence>eps & i<emmax){
    if(link=="logit"){uncure=as.vector(exp(b%*%t(Z))/(1+exp(b%*%t(Z))))}
    if(link=="probit"){uncure=as.vector(pnorm(b%*%t(Z)))}
    if(link=="cloglog"){uncure=as.vector(1-exp(-exp(b%*%t(Z))))}


    s<- psurv(Time,dis,dpar)$basesurv
    h<- psurv(Time,dis,dpar)$basehaz
    survival=drop(s^exp(beta%*%t(X)))
    hhazard=drop(h*exp(beta%*%t(X)))

    #E step

      w=Status*((uncure*mort.s*survival*(mort.h+hhazard))/((1-uncure)*mort.h*mort.s+(uncure*mort.s*survival*(mort.h+hhazard))))+
        (1-Status)*((uncure*mort.s*survival)/((1-uncure)*mort.s+(uncure*mort.s*survival)))


    #M step
    ll1<-function(tryb){
      if(link=="logit"){uncure<- as.vector(exp(tryb%*%t(Z))/(1+exp(tryb%*%t(Z))))}
      if(link=="probit"){uncure<- as.vector(pnorm(tryb%*%t(Z)))}
      if(link=="cloglog"){uncure<- as.vector(1-exp(-exp(tryb%*%t(Z))))}
      -sum(w*log(uncure)+(1-w)*log(1-uncure))
    }
    update_b=try(optim(par=b,fn=ll1))$par


    if(is.matrix(X)){
      qx<- ncol(X)
    }else{qx<- NCOL(X)}

    ll2<-function(par){
      trybeta=par[1:qx]
      if(dis %in% c("exp","weibull","llogis","gamma")){ trydpar<-exp(par[-c(1:qx)])}
      if(dis %in% c("gompertz","lnorm")){ trydpar<-c(par[-c(1:qx)][1], exp(par[-c(1:qx)][-1]))}
      if(dis=="gengamma"){ trydpar<-c(par[-c(1:qx)][1], exp(par[-c(1:qx)][2]), par[-c(1:qx)][3])}
      trys= psurv(Time,dis,trydpar)$basesurv
      tryh= psurv(Time,dis,trydpar)$basehaz
      tryss= drop(trys^exp(trybeta%*%t(X)))
      tryhh= drop(tryh*exp(trybeta%*%t(X)))

      -sum(w*Status*log(mort.h+tryhh)+w*log(tryss))
    }
    update_est=optim(par=c(beta,lgdpar),fn=ll2)$par
    update_beta=update_est[1:qx];update_lgdpar=update_est[-c(1:qx)]

    convergence=sum(c((update_b-b)^2),(c(beta,lgdpar)-update_est)^2)

    b<-update_b
    beta<- update_beta
    lgdpar<- update_lgdpar
    if(dis %in% c("exp","weibull","llogis","gamma")) dpar<-exp(update_lgdpar)
    if(dis %in% c("gompertz","lnorm")) dpar<-c(update_lgdpar[1],exp(update_lgdpar[-1]))
    if(dis=="gengamma") dpar<-c(update_lgdpar[1], exp(update_lgdpar[2]), (update_lgdpar[3]))
    i<- i+1
  }


  #################################
  #variance
  hes.mtx<- diag(rep(0,(length(b)+length(beta)+length(dpar))))
  for(ii in 1:n){
    fst_temp<- grad(com_like,x=c(b,beta,dpar),Time=Time,Status=Status,Z=Z,X=X,dis=dis,
                    w=w,link=link,mort.h=mort.h,mort.s=mort.s,ii=ii)
    tfst_temp<- fst_temp%*%t(fst_temp)
    hes.mtx<- hes.mtx+tfst_temp
  }
  se.est<- sqrt(diag(ginv(hes.mtx)))
  ##################################

  emp <- list(b=b, latencyfit=beta,dpar=dpar,s=s,h=h,Uncureprob=uncure,tau=convergence,se=se.est,w=w)
}

