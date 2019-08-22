#' Background mortality table
#'
#' Background mortality table
#'
#' @param table1 US 2015 Life Table for female
#' @param table2 US 2015 Life Table for male
#' @return None
#'
#'
#' @export
bm_out<- function(table1,table2){
  ageseq=seq(from=0, to=100, by=1)
  trylf<-trylm<-as.list(1:length(ageseq))
  trypf<-trypm<-matrix(rep(0,2*length(ageseq)), ncol=2)
  #lifetable for female with different age
  for(i in 1: length(ageseq)){
  try<-residual_time(table1,table2,sex = "F", age =ageseq[i], year = NULL, country = NULL)
  trylf[[i]]<-try
  trypf[i,]<-exp(optim(c(0,0),lifett,try=try)$par)
  }

  #lifetable for male with different age
  for(i in 1: length(ageseq)){
  try<-residual_time(table1,table2,sex = "M", age =ageseq[i], year = NULL, country = NULL)
  trylm[[i]]<-try
  trypm[i,]<-exp(optim(c(0,0),lifett,try=try)$par)
  }

  try.list<-list(ageseq=ageseq,trypf=trypf,trypm=trypm)
  return(try.list)
}

