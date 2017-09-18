dynfrail_control <- function(opt_fit = TRUE,
                            se = TRUE,
                            se_adj = TRUE,
                            ca_test = TRUE,
                            only_ca_test = FALSE,
                            lik_ci = TRUE,
                            lik_ci_intervals = list(interval = c(-7, 20),
                                                    interval_stable = c(0, 20)),
                            nlm_control = list(),
                            inner_control = list(eps = 0.0001,
                                                 maxit = Inf,
                                                 fast_fit = TRUE,
                                                 verbose = FALSE,
                                                 lower_tol = 20,
                                                 lik_tol = 1)
) {
  # calculate SE as well

  # Here some checks

  if(isTRUE(only_ca_test) & !isTRUE(ca_test))
    stop("control: if only_ca_test is TRUE then ca_test must be TRUE as well")

  if(isTRUE(lik_ci)) {
    if(is.null(lik_ci_intervals))
      stop("opt_control must be a list which contains a named element interval")
    if(length(lik_ci_intervals$interval) != 2)
      stop("interval must be of length 2")
    if(lik_ci_intervals$interval[1] < -7 | lik_ci_intervals$interval[2] > 20)
      warning("extreme values for interval, there might be some numerical trouble")
  }

  inner_c <- function(eps = 0.0001,
                      maxit = Inf,
                      fast_fit = TRUE,
                      verbose = TRUE,
                      lower_tol = 20,
                      lik_tol = 1) {
    list(eps = eps,
         maxit = maxit,
         fast_fit = fast_fit,
         verbose = verbose,
         lower_tol = lower_tol,
         lik_tol = lik_tol)
  }

  inner_control <- do.call(inner_c, inner_control)

  res <- list(opt_fit = opt_fit,
              se = se,
              se_adj = se_adj,
              ca_test = ca_test,
              only_ca_test = only_ca_test,
              lik_ci = lik_ci,
              lik_ci_intervals = lik_ci_intervals,
              nlm_control = nlm_control,
              inner_control = inner_control)
  attr(res, "class") <- c("dynfrail_control")
  res
}




#' Distribution parameters for dynfrail
#'
#' @param dist One of 'gamma', 'stable' or 'pvf'.
#' @param theta A starting value for the 'outer' maximization with respect to the frailty parameter \eqn{\theta}. Must be >0.
#' @param pvfm Only relevant if \code{dist = 'pvf'} is used. It determines which PVF distribution should be used. Must be  larger than -1 and not equal to 0.
#' @param left_truncation Logical. Whether the data set represents left truncated survival times.
#'
#' @return An object of the type \code{emfrail_distribution}, which is mostly used to denote the
#' supported frailty distributions in a consistent way.
#' @export
#'
#' @details The \code{theta} argument must be positive. In the case of gamma or PVF, this is the inverse of
#'  the frailty variance, i.e. the larger the \code{theta} is,
#'  the closer the model is to a Cox model. For the positive stable distribution, the \eqn{\gamma} parameter of the Laplace trnasform is
#'  \eqn{\theta / (1 + \theta)}, with the \eqn{alpha} parameter fixed to 1.
#'
#' @seealso \code{\link{emfrail}, \link{emfrail_control}}
#' @examples
#' emfrail_distribution()
#' # Compound Poisson distribution:
#' emfrail_distribution(dist = 'pvf', theta = 1.5, pvfm = 0.5)
#' # Inverse Gaussian distribution:
#' emfrail_distribution(dist = 'pvf')
dynfrail_distribution <- function(dist = "gamma",
                                  theta = 2,
                                  pvfm = -1/2,
                                  left_truncation = FALSE,
                                  n_ints = NULL,
                                  times = NULL) {

  if (!(dist %in% c("gamma", "stable", "pvf")))
    stop("frailty distribution must be one of gamma, stable, pvf")
  if (length(theta) != 1)
    stop("specify exactly 1 parameter (theta>0) for the frailty")
  if (theta <= 0)
    stop("frailty parameter (theta) must be positive")
  if (dist == "pvf" & (pvfm < -1 | pvfm == 0))
    stop("pvfm must be >-1 and not equal to 0")

  if(!is.logical(left_truncation)) stop("left_truncation must be TRUE or FALSE")
  res <- list(dist = dist, theta = theta, pvfm = pvfm, left_truncation = left_truncation, n_ints = n_ints, times = times)
  attr(res, "class") <- c("dynfrail_distribution")
  return(res)
}



