#' Semiparametric baseline estimation in EM algorithm.
#'
#' Semiparametric baseline estimation in EM algorithm.
#'
#' @param Time observed survival time
#' @param Status right censored status (1=death, 0=censored)
#' @param age the variable of age
#' @param X covariates in latency survival function
#' @param beta parameters for X
#' @param w unobserved cure indicator
#' @param mort.s background survival rate
#' @param mort.h background hazard rate
#' @return None
#'
#'
#'
surv_func<- function(Time,Status,age,X,beta,w,mort.s,mort.h){
  id<- seq(1,length(Time),1)
  temp.data<- data.frame(cbind(id,Time,Status,age,X,w,mort.s,mort.h))
  sort.temp.data<- temp.data[order(temp.data$age),]


  seq_id<- length(id)
  divisor<- 4
  while(divisor<= seq_id){
    if(length(id)%%divisor==0){
      break
    }
    divisor<- divisor+1
  }


  if(divisor<15){
    sdata<- split(sort.temp.data,rep(1:(length(id)/divisor),each=divisor))
  }else{
    while(seq_id>=4){
      if(seq_id %% 4==0){
        break
      }
      seq_id<- seq_id-1
    }

    id_diff<- length(id)-seq_id
    divisor<- 4
    while(divisor<= seq_id){
      if(seq_id%%divisor==0){
        break
      }
      divisor<- divisor+1
    }

    sdata<- split(sort.temp.data[(1:seq_id),],rep(1:(seq_id/divisor),each=divisor))
    sdata[[(seq_id/divisor)+1]]<- sort.temp.data[((seq_id+1):length(id)),]
  }


  base<- list()
  subdata<-list()
  for (i in 1:length(sdata)){
    if(is.matrix(X)){
      XX= as.matrix(sdata[[i]][,5:(5+ncol(X)-1)])
    }else{XX= sdata[[i]]$X}


    base[[i]]<- split_base(Time=sdata[[i]]$Time,Status=sdata[[i]]$Status,
                           X=XX,beta=beta,w=sdata[[i]]$w,
                           mort.s=sdata[[i]]$mort.s,mort.h=sdata[[i]]$mort.h)

    subdata[[i]]<- cbind(sdata[[i]]$id,sdata[[i]]$Time,base[[i]]$basesurv,base[[i]]$tbasehaz)
  }

  comdata<- do.call(rbind,subdata)
  colnames(comdata)<- c("id","Time","basesurv","thaz")
  comdata<- data.frame(comdata)
  final.data<- comdata[order(comdata$id),]
  list(basesurv=final.data$basesurv,haz=final.data$thaz)
}
