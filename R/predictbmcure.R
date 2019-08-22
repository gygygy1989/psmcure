#' Predict survival probabilities.
#'
#' Predict survival probabilities.
#'
#' @param object An object from the function predictbmcure.
#' @param newX The covariates in cure link function.
#' @param newZ The covaraites in latency survival function.
#' @param type The data type can be "whole" for whole data set, "female" for subgroup of female or "male" for subgroup of male.
#' @param \dots Other arguments.
#' @return None
#'
#'
#' @export
predictbmcure <-
  function(object, newX, newZ,type, ...)
  {
    call <- match.call()
    if(!inherits(object, "bmcure")) stop("Object must be bmcure estimate")

    px<-length(object$beta)
    pz<-length(object$b)

    if(is.vector(newZ)) {newZ=as.matrix(newZ)}
    newZ=cbind(1,newZ)
    if(is.vector(newX)) {newX=as.matrix(newX)}

    s0<- drop(split_base(Time=object$Time,Status=object$Status,X=object$X,beta=object$beta,w=object$w,
                    mort.s=object$mort.s,mort.h=object$mort.h)$tsurv) #total baseline survival
    n=length(s0)
    uncure=exp(object$b%*%t(newZ))/(1+exp(object$b%*%t(newZ)))
    cure=1-uncure

    if(type=="whole"){BM<- object$mmort}
    if(type=="male"){BM<- object$m.mmort}
    if(type=="female"){
      BM<- object$f.mmort
      } else {BM<- object$mmort}


    scure=array(0,dim=c(n,nrow(newX)))
    t=array(0,dim=c(n,nrow(newX)))
    spop=array(0,dim=c(n,nrow(newX)))


        ebetaX=exp(object$beta%*%t(newX))
        for(i in 1:nrow(newX)){
        scure[,i]=s0^ebetaX[i]
        }

        for (i in 1:n){
          for (j in 1:nrow(newX)){
        spop[i,j]=uncure[j]*scure[i,j]+(1-uncure[j])*BM
          }
        }


    ttime<- sort(unique(object$Time))
    order.spop= array(0,dim=c(length(ttime),nrow(newX)))

    for(i in 1:length(ttime)){
      for(j in 1:nrow(newX)){
      ppos<- which(object$Time==ttime[i])
      order.spop[i,j]<- max(spop[ppos,j])
      }
    }


    prd<- cbind(order.spop,ttime)
    structure(list(call=call,newcurerate=cure,prediction=prd),class="predictbmcure")
  }
