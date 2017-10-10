#' @export
#' @importFrom stats pchisq printCoefmat
#' @keywords internal
print.dynfrail <- function(x, ...) {

  cat("Call: \n")
  dput(attr(x, "call"))
  cat("\n")

  cat("log-likelihood:", x$loglik[2], "\n")
  cat("theta:", exp(x$logtheta), "|| 1/theta:",  1/exp(x$logtheta), "\n")
  cat("lambda:", exp(x$loglambda), "\n")
  cat("\n")



  if(length(x$coef) > 0) {
    coefmat <- list(
      coef = x$coef,
      "exp(coef)" = exp(x$coef),
      "se(coef)" = sqrt(diag(solve(x$imat))[seq_along(x$coef)]))

    coefmat$z <- coefmat$coef / coefmat$`se(coef)`
    coefmat$p <-  1 - pchisq(coefmat$z^2, df = 1)

    coefmat <- do.call(cbind, coefmat)

    printCoefmat(coefmat)
  }

  invisible(x)

}




