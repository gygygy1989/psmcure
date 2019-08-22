#' Summary of Survival Probabilities from Life Table
#'
#' Summary of Survival Probabilities from Life Table
#'
#' @param table1 US 2015 Life Table for female
#' @param table2 US 2015 Life Table for male
#' @param sex sex
#' @param age age
#' @param year the year of life table
#' @param country the countro of life table
#' @return None
#'
#'
#'
residual_time <- function(table1,table2,sex, age, year, country){
  ## loop over the countries in the trial
  xty <-seq(0,60, by = 0.2)
  Nxty <- length(xty)

  if (sex == "M"){
    scvy_usa <- bm_table( age + xty, table2)
  }else{
    scvy_usa <- bm_table( age + xty, table1)
  }
  ret_obj <- list()
  #mean survival curve
  ret_obj$surv_mean_usa <- scvy_usa
  ret_obj$time <- xty
  return(ret_obj)
}

