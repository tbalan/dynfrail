#' Control parameters for dynfrail
#'
#' @param opt_fit Logical. Whether the outer optimization should be carried out.
#' If \code{FALSE}, then the frailty parameter is treated as fixed and the \code{emfrail} function returns only log-likelihood. See details.
#' @param se Ignored
#' @param se_adj Ignored
#' @param nlm_control A list of named arguments to be sent to \code{nlm} for the outer optimization.
#' @param inner_control A list of parameters for the inner optimization. See details.
#'
#' @return An object of the type \code{emfrail_control}.
#' @export
#'
#' @details
#' The \code{nlm_control} argument should not overalp with \code{hessian}, \code{f} or \code{p}.
#'
#' The \code{inner_control} argument should be a list with the following items:
#' \itemize{
#' \item{\code{eps}}{ A criterion for convergence of the EM algorithm (difference between two consecutive values of the log-likelihood)}
#' \item{\code{maxit}}{ The maximum number of iterations between the E step and the M step}
#' \item{\code{verbose}}{ Logical, whether details of the optimization should be printed}
#' \item{\code{lik_tol}}{ For values higher than this, the algorithm returns a warning when the log-likelihood decreases between EM steps. Technically, this should not happen, but
#' if the parameter \eqn{\theta} is somewhere really far from the maximum, numerical problems might lead in very small likelihood decreases.
#' }}
#'
#' The starting value of the outer optimization may be set in the \code{dynfrail_dist()} argument.
#'
#' @seealso \code{\link{emfrail}}, \code{\link{emfrail_dist}}, \code{\link{emfrail_pll}}
#' @examples
#' dynfrail_control()
#' dynfrail_control(inner_control = list(maxit = 2))
#'
dynfrail_control <- function(opt_fit = TRUE,
                            se = TRUE,
                            se_adj = TRUE,
                            nlm_control = list(),
                            inner_control = list(eps = 0.0001,
                                                 maxit = 100,
                                                 verbose = TRUE,
                                                 # lower_tol = 20,
                                                 lik_tol = 1)
) {
  # calculate SE as well


  inner_c <- function(eps = 0.0001,
                      maxit = Inf,
                      verbose = TRUE,
                      lik_tol = 1) {
    list(eps = eps,
         maxit = maxit,
         verbose = verbose,
         lik_tol = lik_tol)
    }

  inner_control <- do.call(inner_c, inner_control)

  res <- list(opt_fit = opt_fit,
              se = se,
              se_adj = se_adj,
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
                                  lambda = 0.1,
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

  # if(!is.logical(left_truncation)) stop("left_truncation must be TRUE or FALSE")
  res <- list(dist = dist, theta = theta, pvfm = pvfm, lambda = lambda, n_ints = n_ints, times = times)
  attr(res, "class") <- c("dynfrail_distribution")
  return(res)
}



