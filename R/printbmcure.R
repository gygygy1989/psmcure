#' Print out estimation.
#'
#' Print out estimation.
#'
#' @param x An object from estimation in the function bmcure.
#' @param \dots Other plot arguments.
#' @return None
#'
#'
#' @export
printbmcure <-
  function(x, ...)
  {
    if(!is.null(cl <- x$call)) {
      cat("Call:\n")
      dput(cl)
    }
    cat("Cure Estimate:\n")
      b <- array(x$b,c(length(x$b),4))
      rownames(b) <- x$bnm
      colnames(b) <- c("Estimate","Std.Error","Z value","Pr(>|Z|)")
      b[,2] <- x$bse
      b[,3] <- x$bzvalue
      b[,4] <- x$bpvalue


    print(b)
    cat(paste0("Arm 1 cure rate=",round(x$cure[1]*100,3),"%"," ; ","Arm 2 cure rate=",round(x$cure[2]*100,3),"%","\n"))
    cat("\n")

    cat("Latency Survival Model:\n")
      beta <- array(x$beta,c(length(x$beta),4))
      rownames(beta) <- x$betanm
      colnames(beta) <- c("Estimate","Std.Error","Z value","Pr(>|Z|)")
      beta[,2] <- x$betase
      beta[,3] <- x$betazvalue
      beta[,4] <- x$betapvalue
    print(beta)
    cat("-----------------------------------------------------------------\n")
    invisible(x)
  }
