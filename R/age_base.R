#' Baseline estimation stratified by age
#'
#' Baseline estimation stratified by age
#'
#' @param Time observed survival time
#' @param Status right censored status (1=death, 0=censored)
#' @param X covariates in latency survival function
#' @param beta parameters for X
#' @param w unobserved cure indicator
#' @param mort.s background survival rate
#' @param mort.h background hazard rate
#' @return None
#'
#'
#'
split_base<- function(Time,Status,X,beta,w,mort.s,mort.h){

  nn<- length(Time)
  death_point<-sort(unique(subset(Time,Status==1)))
  lambda<-numeric()
  event<-numeric()

  if(length(death_point)==0){
    #total baseline haz
    tbasehaz=rep(0,nn)

    #total cumulative baseline haz
    tHhaz=rep(0,nn)

    #total baseline surv
    tsurv=rep(1,nn)

    #baseline surv
    basesurv<- tsurv/exp(log(mort.s)/drop(exp(beta%*%t(X))))
  }else{

    for (i in 1:length(death_point)){
      event[i]<-sum(Status*(w)*as.numeric(Time==death_point[i]))
      temp=sum(as.numeric(Time>=death_point[i])*(w)*drop(exp(beta%*%t(X))))
      temp1<-event[i]
      lambda[i]<-temp1/temp
    }


    tbasehaz<- numeric()
    for (i in 1:length(Time)){
      tbasehaz[i]<- sum(as.numeric(Time[i]==death_point)*lambda)
    }


    tHhaz<-numeric()
    for(i in 1:length(Time)){
      tHhaz[i]<- sum(as.numeric(Time[i]>=death_point)*lambda)
      if (Time[i]>max(death_point)) tHhaz[i]<- Inf  #zerotail constrain
      if (Time[i]<min(death_point)) tHhaz[i]<- 0
    }

    tsurv<- exp(-tHhaz)

    basesurv<- tsurv/exp(log(mort.s)/drop(exp(beta%*%t(X))))
  }

  list(tsurv=tsurv,basesurv=basesurv,tbasehaz=tbasehaz)
}
