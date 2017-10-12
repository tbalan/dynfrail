#' Control parameters for dynfrail
#'
#' @param nlm_control A list of named arguments to be sent to \code{nlm} for the outer optimization.
#' @param inner_control A list of parameters for the inner optimization. See details.
#'
#' @return An object of the type \code{dynfrail_control}.
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
#' @seealso \code{\link{dynfrail}}, \code{\link{dynfrail_dist}}
#' @examples
#' dynfrail_control()
#' # this stops each EM (inner maximization) after 10 iterations, event if it did not
#' # reach the maximum.
#' dynfrail_control(inner_control = list(maxit = 10))
#'
dynfrail_control <- function(nlm_control = list(stepmax = 1),
                            inner_control = list(eps = 0.0001,
                                                 maxit = 100,
                                                 verbose = FALSE,
                                                 # lower_tol = 20,
                                                 lik_tol = 1)
) {
  # calculate SE as well


  inner_c <- function(eps = 0.0001,
                      maxit = Inf,
                      verbose = FALSE,
                      lik_tol = 1) {
    list(eps = eps,
         maxit = maxit,
         verbose = verbose,
         lik_tol = lik_tol)
    }

  inner_control <- do.call(inner_c, inner_control)

  res <- list(nlm_control = nlm_control,
              inner_control = inner_control)
  attr(res, "class") <- c("dynfrail_control")
  res
}




#' Distribution parameters for dynfrail
#'
#' @param dist One of 'gamma', 'stable' or 'pvf'.
#' @param theta Frailty distribution parameter. Must be >0.
#' @param pvfm Only relevant if \code{dist = 'pvf'} is used. It determines which PVF distribution should be used. Must be larger than -1 and not equal to 0.
#' @param lambda Frailty autocorrelation parameter. Must be >0.
#' @param n_ints For piece-wise constant frailty, the number of intervals. With \code{n_ints = 0}, the classical shared frailty scenario is obtained.
#' @param times A vector of time points which determine the piecewise-constant interval for the frailty. Overrides \code{n_ints}.
#' @return An object of the type \code{dynfrail_dist}, which is mostly used to denote the
#' supported frailty distributions in a consistent way.
#' @export
#'
#' @details The \code{theta} and \code{lambda} arguments must be positive. In the case of gamma or PVF, \code{theta} is the inverse of
#'  the frailty variance, i.e. the larger the \code{theta} is,
#'  the closer the model is to a Cox model. When \code{dist = "pvf"} and \code{pvfm = -0.5}, the inverse Gaussian
#'  distribution is obtained.
#'  For the positive stable distribution, the \eqn{\gamma} parameter of the Laplace transform is
#'  \eqn{\theta / (1 + \theta)}, with the \eqn{alpha} parameter fixed to 1.
#'
#' @seealso \code{\link{dynfrail}, \link{dynfrail_control}}
#' @examples
#' dynfrail_dist()
#' # Compound Poisson distribution:
#' dynfrail_dist(dist = 'pvf', theta = 1.5, pvfm = 0.5)
#' # Inverse Gaussian distribution:
#' dynfrail_dist(dist = 'pvf')
dynfrail_dist <- function(dist = "gamma",
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
  attr(res, "class") <- c("dynfrail_dist")
  return(res)
}



