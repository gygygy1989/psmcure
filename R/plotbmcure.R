#' Plot the predicted survival function.
#'
#' Plot the predicted survival function.
#'
#' @param object An object from the function predictbmcure.
#' @param type The type of plot, such as "S" or "l".
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @param \dots Other plot arguments.
#' @return None
#'
#'
#' @export
plotbmcure <-
  function(object, type="S", xlab="Time",ylab="Predicted Survival Probability", ...)
  {
    predsurv <- object$prediction

      if(length(object$newuncurerate)==1) {
        plot(predsurv[,ncol(predsurv)],predsurv[,1], type="l",ylim=c(0,1))
      } else {
      matplot(predsurv[,ncol(predsurv)],predsurv[,1:(ncol(predsurv)-1)],
                col=1,type="l",lty=1:(ncol(predsurv)-1),lwd=1,
                xlab=xlab,ylab=ylab,ylim=c(0,1))
      }
  }
