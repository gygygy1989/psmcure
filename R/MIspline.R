#' Basis Splines in Flexible Parametric Method.
#'
#' Basis Splines in Flexible Parametric Method.
#'
#' @param x a row vector (grid of survival time)
#' @param order degree of splines
#' @param knots a sequence of increasing points
#' @return None
#'
#'
#'
MIspline<-function(x,order,knots){
  # get M spline and I spline matrix with order
  # x is a row vector
  # k is the order of I spline
  # knots are a sequence of increasing points
  # the number of free parameters in M spline is the length of knots plus 1.
  # Created by Bo Cai and Lianming Wang, Oct. 2009, Sept. 2011

  ### get Mspline bases ###
  k1=order
  m=length(knots)
  n1=m-2+k1 # number of parameters
  t1=c(rep(1,k1)*knots[1], knots[2:(m-1)], rep(1,k1)*knots[m]) # newknots

  tem1=array(rep(0,(n1+k1-1)*length(x)),dim=c(n1+k1-1, length(x)))
  for (l in k1:n1){
    tem1[l,]=(x>=t1[l] & x<t1[l+1])/(t1[l+1]-t1[l])
  }

  if (order==1){
    mbases=tem1
  }else{
    mbases=tem1
    for (ii in 1:(order-1)){
      tem=array(rep(0,(n1+k1-1-ii)*length(x)),dim=c(n1+k1-1-ii, length(x)))
      for (i in (k1-ii):n1){
        tem[i,]=(ii+1)*((x-t1[i])*mbases[i,]+(t1[i+ii+1]-x)*mbases[i+1,])/(t1[i+ii+1]-t1[i])/ii
      }
      mbases=tem
    }
  }

  ### get Ispline bases ###
  k=order+1
  n=m-2+k # number of parameters
  t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots

  yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
  for (l in k:n){
    yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
  }

  yytem1=yy1
  for (ii in 1:order){
    yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
    for (i in (k-ii):n){
      yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
    }
    yytem1=yytem2
  }


  index=rep(0,length(x))
  for (i in 1:length(x)){
    index[i]=sum(t<=x[i])
  }

  ibases=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))

  if (order==1){
    for (i in 2:n){
      ibases[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
    }
  }else{
    for (j in 1:length(x)){
      for (i in 2:n){
        if (i<(index[j]-order+1)){
          ibases[i-1,j]=1
        }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
          ibases[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
        }else{
          ibases[i-1,j]=0
        }
      }
    }
  }
  return(list(mbases,ibases))
}
