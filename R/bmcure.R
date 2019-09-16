#' Main function to estimate parameters and print out
#'
#' Main function to estimate parameters and print out
#'
#' @param formula A survival formula based on the Surv() function, including survival time, right censoring indicator and covariates.
#' @param cureform The covariates included in the cure link function.
#' @param link The type of cure rate link function, including "logit","probit", and "cloglog".
#' @param data The data set needed.
#' @param table1 US 2015 Life Table for female.
#' @param table2 US 2015 Life Table for male.
#' @param na.action The action to deal with NA/missing values.
#' @param method The method to estimate parameters, including "semiparametric","flexible" and "parametric".
#' @param dis If "parametric" is specified for method, then dis must be provided, such as "exp","weibull","llogis","gamma","gompertz","lnorm" and "gengamma".
#' @importFrom stats na.omit model.frame model.matrix model.extract pweibull dweibull
#' @return None
#'
#'
#' @export
bmcure<- function(formula,cureform,link,data,table1,table2,na.action=na.omit,method="parametric",dis="weibull"){
  call<- match.call()
  cat("running... please wait \n")
  data<- na.action(data)
  n<- nrow(data)
  mf<- model.frame(formula,data)
  curevar<- all.vars(cureform)

  #########################
  #match life table
  try.list<- bm_out(table1,table2)
  trypf=try.list$trypf;trypm=try.list$trypm;ageseq=try.list$ageseq
  nn<- length(ageseq)
  sex1=rep(1,nn);sex0=rep(0,nn) #1 is female;0 is male
  fdata<-data.frame(cbind(ageseq,sex1,trypf)); mdata<-data.frame(cbind(ageseq,sex0,trypm))
  data$id=seq(1,length(data$age),1)
  n= length(data$age)

  par1<- vector()
  par2<- vector()
  for (i in 1:n){
    if (data$sex[i]==2){ #2 is female; 1 is male
      pos=which(fdata$ageseq==data$age[i])
      par1[i]=fdata$V3[pos]
      par2[i]=fdata$V4[pos]
    }
    if(data$sex[i]==1){
      pos=which(mdata$ageseq==data$age[i])
      par1[i]=mdata$V3[pos]
      par2[i]=mdata$V4[pos]
    }
  }

  data$mort.s<- 1-pweibull(data$Time+1e-10,par1,par2)+1e-10
  data$mort.h<- dweibull(data$Time+1e-10,par1,par2)/(data$mort.s)
  mmort<- mean(data$mort.s)
  m.data<- data[which(data$sex==1),]
  f.data<- data[which(data$sex==2),] #2 is female
  m.mmort<- mean(m.data$mort.s)
  f.mmort<- mean(f.data$mort.s)
  #########################################

  Z<- as.matrix(cbind(rep(1,n),data[,curevar]))
  colnames(Z)<- c("Intercept",curevar)

  Y<- model.extract(mf,"response")
  X<- model.matrix(attr(mf,"terms"),mf)
  if(!inherits(Y,"Surv")) stop("Response must be time")
  Time<- Y[,1]
  Status<- Y[,2]
  age<- data[,"age"]
  bnm<- colnames(Z)
  nb<- ncol(Z)

    betanm<- colnames(X)[-1]
    nbeta<- ncol(X)-1

  mort.s<- data$mort.s
  mort.h<- data$mort.h
  X<- X[,-1]

  if (method=="semiparametric"){
    emfit<- emback(Time,Status,age,Z,X,mort.s,mort.h,link=link,emmax=100,eps=1e-7)
    se<- emfit$se
  }
  if(method=="flexible"){
    # initial values
    w= Status
    n= length(Status)
    # initial b
    b=glm(w~Z[,-1],family=quasibinomial(link=link))$coef
    # initial beta
    beta=coxph(Surv(Time,Status)~X+offset(log(w+1e-10)),method="breslow")$coef
    ##########################
    # select knots
    grids=Time
    knots<- knot.fun(Time,Status,Z,X,b,beta,w,mort.s,mort.h,link=link,emmax=100,eps=1e-7,maxk=10)
    MI=MIspline(grids,order=3,knots) #order=2 or 3
    Ibigs=MI[[2]]   # I spline basis evaluated at grids
    Mbigs=MI[[1]]  # M spline basis evaluated at grids
    mm<- nrow(Ibigs)
    r<- rep(0.1,mm)

    emfit<- flex.spline(Time,Status,Z,X,Mbigs,Ibigs,r,b,beta,w,mort.s,mort.h,link="logit",emmax=100,eps=1e-10)
    b<- emfit$b
    beta<- emfit$latency
    r<- emfit$r
    w<- emfit$w
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
    se<- se0r[1:(length(b)+length(beta))]
  }
  if (method=="parametric"){
    emfit<- emp(Time,Status,age,Z,X,mort.s,mort.h,link="logit",emmax=100,eps=1e-7,dis=dis)
    se<- emfit$se
  }
  b<- emfit$b
  beta<- emfit$latencyfit
  bse<- se[1:nb]
  betase<- se[(nb+1):(nb+nbeta)]
  s<- emfit$s
  h<- emfit$h
  cure<- unique(1-emfit$Uncureprob)
  w<- emfit$w

  fit<- list()
  class(fit)<- c("bmcure")
  fit$b<- b
  fit$cure<- cure
  fit$beta<- beta
  fit$bse<- bse
  fit$bzvalue<- b/bse
  fit$bpvalue<- (1-pnorm(abs(fit$bzvalue)))*2
  fit$betase<- betase
  fit$betazvalue<- beta/betase
  fit$betapvalue<- (1-pnorm(abs(fit$betazvalue)))*2

  #following is needed to predict
  fit$Time<- Time
  fit$Status<- Status
  fit$X<- X
  fit$w<- w
  fit$mort.s<- mort.s
  fit$mort.h<- mort.h
  fit$mmort<- mmort
  fit$m.mmort<- m.mmort
  fit$f.mmort<- f.mmort
  cat("done.\n")
  fit$call<- call

  printbmcure(fit)

}
