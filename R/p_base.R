#' Parametric baseline estimation.
#'
#' Parametric baseline estimation.
#'
#' @param Time Survival time from data set.
#' @param dis The parametric distribution such as "exp","weibull","llogis","gamma","gompertz","lnorm" and "gengamma".
#' @param dpar the parameters for different distributions.
#' @importFrom flexsurv hexp hweibull pllogis hllogis hgamma pgompertz hgompertz pgengamma hgengamma hlnorm
#' @return None
#'
#'
#'
psurv<-function(Time,dis,dpar){
  if(dis=="exp"){
    s<- 1-pexp(Time,dpar[1])
    h<- hexp(Time,dpar[1])
  }
  if(dis=="weibull"){
    s<- 1-pweibull(Time,dpar[1],dpar[2])
    h<- hweibull(Time,dpar[1],dpar[2])
  }
  if(dis=="llogis"){
    s<- 1-pllogis(Time,dpar[1],dpar[2])
    h<- hllogis(Time,dpar[1],dpar[2])
  }
  if(dis=="gamma"){
    s<- 1-pgamma(Time,dpar[1],dpar[2])
    h<- hgamma(Time,dpar[1],dpar[2])
  }
  if(dis=="gompertz"){
    s<- 1-pgompertz(Time, dpar[1], dpar[2])
    h<- hgompertz(Time, dpar[1], dpar[2])
  }
  if(dis=="gengamma"){
    s<- 1-pgengamma(Time, dpar[1], dpar[2], dpar[3])
    h<- hgengamma(Time, dpar[1], dpar[2], dpar[3])
  }

  if(dis=="lnorm"){
    s<- 1-plnorm(Time, dpar[1], dpar[2])
    h<- hlnorm(Time, dpar[1], dpar[2])
  }


  list(basesurv=s,basehaz=h)
}
