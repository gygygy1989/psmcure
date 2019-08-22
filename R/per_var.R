#' Perturbation variance estimation
#'
#' Perturbation variance estimation
#'
#' @param b parameters in cure link function.
#' @param Z The covariates included in the cure link function.
#' @param link The type of cure rate link function, including "logit","probit", and "cloglog".
#' @param Status Right censoring status: 1 is death, 0 is censored.
#' @param w unobserved cure status indicator.
#' @param beta The estimated parameters in latency survival function.
#' @param s The estimated baseline survival rate.
#' @param th  The estiamted total baseline hazard rate.
#' @param X The covaraites for latency survival function.
#' @param mort.h background hazard rate from life table.
#' @param mort.s background survival rate from life table.
#'
#' @importFrom Matrix nearPD
#' @importFrom numDeriv grad
#'
#' @return None
#'
#'
#' @export
per_fun<- function(b=b,Z=Z,link=link,Status=Status,w=w,beta=beta,s=s,th=th,X=X,mort.h=mort.h,mort.s=mort.s){
  n<- length(Status)
  part1<- matrix(nrow=n,ncol=length(b))
  part2<- matrix(nrow=n,ncol=length(beta))

  for(ii in 1:n){
    part1[ii,]<- grad(fish_func1,x=b,Z=Z,w=w,link=link,ii=ii)
  }

  d1=0.1
  ind<- diag(1,nrow=length(beta),ncol=length(beta))

  for(ii in 1:n){
    for(jj in 1:length(beta)){
      ll1<- fish_func2(beta=(beta-2*d1*ind[jj,]),s=s,th=th,Status=Status,X=X,w=w,mort.h=mort.h,mort.s=mort.s,ii=ii)
      ll2<- fish_func2(beta=(beta-d1*ind[jj,]),s=s,th=th,Status=Status,X=X,w=w,mort.h=mort.h,mort.s=mort.s,ii=ii)
      ll3<- fish_func2(beta=(beta+d1*ind[jj,]),s=s,th=th,Status=Status,X=X,w=w,mort.h=mort.h,mort.s=mort.s,ii=ii)
      ll4<- fish_func2(beta=(beta+2*d1*ind[jj,]),s=s,th=th,Status=Status,X=X,w=w,mort.h=mort.h,mort.s=mort.s,ii=ii)
      part2[ii,jj]<- (1/(12*d1))*(ll1-8*ll2+8*ll3-ll4)
    }
  }

  s.func<-matrix(nrow=length(Status),ncol=length(c(b,beta)))
  i.func<-diag(rep(0,length(c(b,beta))))
  for (ii in 1:length(Status)){
    s.func[ii,]<- cbind(part1,part2)[ii,]
    temp.func<- s.func[ii,]%*%t(s.func[ii,])
    i.func<- i.func+temp.func
  }

  se.est<- sqrt(diag(nearPD(ginv(i.func))$mat))
  return(se.est)
}





