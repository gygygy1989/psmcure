#' Knots Selection in Flexible Parametric EM Algorithm.
#'
#' Knots Selection in Flexible Parametric EM Algorithm.
#'
#' @param Time Survival time from data set.
#' @param Status Right censoring status: 1 is death, 0 is censored.
#' @param Z The covariates for cure link function.
#' @param X The covaraites for latency survival function.
#' @param b Initial value for b
#' @param beta  Initial value for beta
#' @param w Initial value for w
#' @param mort.s background survival rate from life table.
#' @param mort.h background hazard rate from life table.
#' @param link The type of cure rate link function, including "logit","probit", and "cloglog".
#' @param emmax The maximum iteration is 100.
#' @param eps The stopping criteria for convergence is 1e-7.
#' @param maxk The maximum number of knots to search
#' @importFrom stats glm coef pnorm optim quantile
#' @importFrom numDeriv grad
#' @importFrom MASS ginv
#' @return None
#'
#'
#'
#'
knot.fun<- function(Time,Status,Z,X,b,beta,w,mort.s,mort.h,link="logit",emmax,eps,maxk){ #select knots
  grids=Time
  aic.out<- c()
  count<- c()
  ti<- Time[Status==1]
  min_t<- min(ti)+(1e-7)
  max_t<- max(ti)+(1e-7)

  for(k in 3:maxk){ #i is number of knots
    id<- seq(0,1,length.out=k)
    id<- id[-c(1,k)]
    knots<- c(min_t,quantile(ti,id),max_t)

    MI=MIspline(grids,order=3,knots) #order=2 or 3
    Ibigs=MI[[2]]   # I spline basis evaluated at grids
    Mbigs=MI[[1]]  # M spline basis evaluated at grids
    mm<- nrow(Ibigs)
    r<- rep(0.1,mm)
    out<- flex.spline(Time,Status,Z,X,Mbigs,Ibigs,r,b,beta,w,mort.s,mort.h,link="logit",emmax=100,eps=1e-7)
    out_est<- c(out$b,out$latencyfit,out$r)
    logl<- like.com(est0r=out_est,link,X,Z,Mbigs,Ibigs,Status,w=out$w,mort.s,mort.h)
    aic<- -2*logl+2*length(out_est)
    aic.out<- c(aic.out,aic)
    count<- c(count,k)
  }
  pos<- which(aic.out==min(aic.out))
  k<- seq(3,maxk)
  id_pos<- seq(0,1,length.out=k[pos])
  id_pos<- id_pos[-c(1,k[pos])]
  knots<- c(min_t,quantile(ti,id_pos),max_t)
  #cat("#knots=",k[pos],"\n")
  return(knots)
}
