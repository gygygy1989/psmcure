#' Extract Survival Probabilities from Life Table
#'
#' Extract Survival Probabilities from Life Table
#'
#' @param years age
#' @param table life table
#' @return None
#'
#'
#'
bm_table <- function(years,table){
  n <- length(years);
  v_ret <- rep(0, length(years))
  for (i in 1:n){
    yr_l <- floor(years[i])
    yr_h <- ceiling(years[i])
    if ( yr_l != yr_h){
      s_h <- table$lx[yr_h+1]/1e5
      s_l <- table$lx[yr_h+0]/1e5 #number surviving to age x
      value <- (s_h - s_l)*(years[i] - yr_l) + s_l
    }else{
      value <- table$lx[yr_h+1]/1e5
    }
    if (yr_h > 100){ value <- 0 }
    v_ret[i]<- value
  }
  return(v_ret)
}
