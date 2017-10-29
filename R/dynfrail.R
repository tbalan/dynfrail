#' Fitting dynamic frailty models with the EM algorithm
#'
#'
#' @importFrom survival Surv survSplit
#' @importFrom magrittr "%>%"
#' @importFrom Rcpp evalCpp
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr separate
#' @importFrom dplyr arrange_
#' @importFrom stats nlm drop.terms model.frame quantile terms model.matrix
#' @useDynLib dynfrail, .registration=TRUE
#' @include dynfrail_aux.R
#' @include dynfrail_fit.R
#'
#' @param formula A formula that contains on the left hand side an object of the type \code{Surv}
#' and on the right hand side a \code{+cluster(id)} statement.
#' @param data A data frame in which the formula argument can be evaluated
#' @param distribution An object as created by \code{\link{dynfrail_dist}}
#' @param control An object as created by \code{\link{dynfrail_control}}
#' @param ... Other arguments, currently used to warn about deprecated argument names
#' @export
#'
#' @return A \code{dynfrail} object that contains the following fields:
#' \item{coefficients}{A named vector of the estimated regression coefficients}
#' \item{hazard}{The breslow estimate of the baseline hazard at each event time point, in chronological order}
#' \item{imat}{Fisher's information matrix corresponding to the coefficients and hazard, assuming \eqn{\theta, \lambda} constant}
#' \item{logtheta}{The point estimate of the logarithm of the frailty parameter \eqn{\theta}. See details.}
#' \item{loglambda}{The point estimate of the logarithm of the autocorrelation parameter \eqn{\lambda}. See details.}
#' \item{frail}{A \code{data.frame} containing the variables: \code{id} (cluster id), \code{interval} (for piecewise constant frailty, the label of the interval
#' on which the frailty is constant), \code{Y} (a \code{Surv} object which determines a starting and a stopping time for each row), \code{frail} (the empirical
#' Bayes estimates of the piecewise constant frailty corresponding to that specific cluster and that specific time period)}
#' \item{tev}{The time points of the events in the data set, this is the same length as hazard}
#' \item{loglik}{A vector of length two with the log-likelihood of the starting Cox model
#' and the maximized log-likelihood}
#' \item{formula}{The original formula argument}
#' \item{distribution}{The original distribution argument}
#' \item{control}{The original control argument}
#' @export
#'
#' @references Putter, H., & Van Houwelingen, H. C. (2015). Dynamic frailty models based on compound birth–death processes. Biostatistics, 16(3), 550-564.
#'
#' @details This function fits dynamic frailty models where the intensity of the process is described by
#' \deqn{\lambda(t) = Z(t) \exp(\beta^\top x) \lambda_0(t).} As in regular frailty models, the random effect is
#' shared by observations from a cluster, or by recurrent event episodes within an individual.
#' This implementation generally follows the lines of Putter & van Houwelingen (2015). The maximum likelihood
#' estimates are obtained with an exact E step.
#'
#' \eqn{Z(t)} has two parameters: \eqn{\theta} plays the role of the spread of the frailty distribution. For the frailty distributons with finite
#' variance (all except the positive stable) this is the inverse of the variance, so that 0 corresponds to infinite variance and infinity to
#' variance 0. The second parameter \eqn{\lambda} determines how much variation in time is in \eqn{Z(t)}, so that
#' \deqn{cor(Z(t_1), Z(t_2)) = exp(-\lambda (t_2 - t_1)).} Note that this heavily depends on the time scale, so
#' the starting value in the \code{distribution} should reflect that.
#'
#' By default, the program must calculate \eqn{Z(t)} for each cluster and for each event time point in the data.
#' This is computationally challenging. An option is to use the \code{nints} argument in the \code{control} argument.
#' This considers \eqn{Z(t)} to be piecewise constant over \code{nints + 1} intervals. These intervals are determined
#' automatically so that there are roughly an equal number of observations for each interval. Using \code{nints = 0}
#' is equivalent to fitting a shared frailty model with the \code{frailtyEM} package.
#'
#' It is recommended that the user starts with \code{nints = 0} and then slowly increase the number of intervals.
#' Other options for performance may be set within the \code{control} argument. Also, this could be tried out first
#' on a subset of the data.
#'
#' For computational reasons, the standard errors of \eqn{\theta} and \eqn{\lambda} are not calculated, and the standard errors of
#' the regression coefficients are obtained under the assumption that the frailty distribution is fixed.
#'
#' @examples
#' # 5 piecewise constant intervals
#' m2 <- dynfrail(Surv(time, status) ~ rx + sex + cluster(litter),
#' data = rats,
#' distribution = dynfrail_dist(n_ints = 4))
#'
#' \dontrun{
#' #' # essentially a gamma frailty fit
#' m1 <- dynfrail(Surv(time, status) ~ rx + sex + cluster(litter),
#' data = rats,
#' distribution = dynfrail_dist(n_ints = 0))
#'
#' # completely semiparametric gamma frailty
#' m2 <- dynfrail(Surv(time, status) ~ rx + sex + cluster(litter),
#' data = rats)
#' }
dynfrail <- function(formula, data,
                     distribution = dynfrail_dist(),
                     control = dynfrail_control(),
                     ...)
{


  if(!inherits(distribution, "dynfrail_dist"))
    stop("distribution argument misspecified; see ?dynfrail_dist()")

  if(!inherits(control, "dynfrail_control"))
    stop("control argument misspecified; see ?dynfrail_control()")


  Call <- match.call()

  # Here checking inputs
  if(missing(formula)  | missing(data)) stop("Missing arguments")

  cluster <- function(x) x

  Y <- model.frame(formula, data)[[1]]

  if(!inherits(Y, "Surv")) stop("left hand side not a survival object")
  if(ncol(Y) != 3) {
    Y <- Surv(rep(0, nrow(Y)), Y[,1], Y[,2])
  }

  tev_unique_ord <- sort(unique(Y[,2][Y[,3]==1]))

  # browser()

  if(!is.null(distribution$n_ints)) {

    quants <- quantile(tev_unique_ord,
             probs = seq(from = 0, to = 1,
                         length.out = distribution$n_ints + 2), type = 1)


    cut <- quants[-c(1, length(quants))]

    # this is a bit dodgy but this is only needed further to
    # assign time points to intervals of the frailty
    tev_unique_ord <- quants[-1]

  } else
    if(!is.null(distribution$times)) cut <- distribution$times else
      cut <- tev_unique_ord[-length(tev_unique_ord)]

  df_dynfrail <- survSplit(formula, data = data, cut = cut, episode = "interval_")


  names(df_dynfrail)[grep("cluster", names(df_dynfrail))] <- "id_"


  pos_id <- grep("cluster", attr(terms(formula), "term.labels"))

  terms2 <- drop.terms(terms(formula), dropx = pos_id, keep.response = TRUE)

  mf <- model.frame(terms2, df_dynfrail)

  Y <- mf[[1]]
  if(ncol(Y) != 3) {
    Y <- Surv(df_dynfrail$tstart, Y[,1], Y[,2])
  }
  # get the model matrix
  X1 <- model.matrix(terms2, df_dynfrail)

  X <- X1[,-c(1), drop=FALSE]
  # note: X has no attributes, in coxph it does.

  mcox <- survival::agreg.fit(x = X,
                              y = Y,
                              strata = NULL, offset = NULL, init = NULL,
                              control = survival::coxph.control(),
                              weights = NULL,
                              method = "breslow",
                              rownames = NULL)

  # Now we need to get the hazard out of this for every row

  if(length(X) == 0) {
    newrisk <- 1
    exp_g_x <- matrix(rep(1, length(mcox$linear.predictors)), nrow = 1)
    g <- 0
    g_x <- t(matrix(rep(0, length(mcox$linear.predictors)), nrow = 1))

  } else {
    x2 <- matrix(rep(0, ncol(X)), nrow = 1, dimnames = list(123, dimnames(X)[[2]]))
    x2 <- (scale(x2, center = mcox$means, scale = FALSE))
    newrisk <- exp(c(x2 %*% mcox$coefficients) + 0)
    exp_g_x <- exp(mcox$coefficients %*% t(X))
    g <- mcox$coefficients
    g_x <- t(mcox$coefficients %*% t(X))

  }


  # what I need here is
  explp <- exp(mcox$linear.predictors) # these are with centered covariates

  #nev_id <- rowsum(Y[,3], id) # nevent per id or am I going crazy

  nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))
  esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))

  # the stuff that won't change
  death <- (Y[, ncol(Y)] == 1)
  nevent <- as.vector(rowsum(1 * death, Y[, ncol(Y) - 1])) # per time point
  time <- sort(unique(Y[,2])) # unique tstops incl. censoring times

  # this gives the next entry time for each unique tstop (not only event)
  etime <- c(0, sort(unique(Y[, 1])),  max(Y[, 1]) + min(diff(time)))
  indx <- findInterval(time, etime, left.open = TRUE) # left.open  = TRUE is very important

  # this gives for every tstart (line variable) after which event time did it come
  # indx2 <- findInterval(Y[,1], time, left.open = FALSE, rightmost.closed = TRUE)
  indx2 <- findInterval(Y[,1], time)

  time_to_stop <- match(Y[,2], time)
  # order_id <- findInterval(id, unique(id))

  atrisk <- list(death = death, nevent = nevent,
                 #nev_id = nev_id,
                 #order_id = order_id,
                 time = time, indx = indx, indx2 = indx2,
                 time_to_stop = time_to_stop)
  atrisk$x2 <- x2

  # Some stuff that won't change
  atrisk$id_interval <- paste0(df_dynfrail$id_, "_",df_dynfrail$interval_)

  # events / interval for each id
  death_id_interval <- rowsum(Y[,3], group = atrisk$id_interval, reorder = TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::separate("rowname", into = c("id", "interval"), sep = "_", convert = TRUE) %>%
    dplyr::arrange_("id", "interval")

  delta <- split(death_id_interval$V1, death_id_interval$id)
  intervals <- split(death_id_interval$interval, death_id_interval$id)

  atrisk$times_incluster <- lapply(intervals, function(x)
    tev_unique_ord[x]
  )

  # for each cluster, in which interval do the events fall
  atrisk$events_incluster <- lapply(delta, function(x) {
    rep(1:length(x), x)
  })

  # for each row, to which frailty interval it belongs
  atrisk$interval_incluster <-
    split(df_dynfrail$interval_, df_dynfrail$id_)

  # First calculation of the cumulative hazard
  nrisk <- nrisk - c(esum, 0,0)[indx]

  haz <- nevent/nrisk * newrisk

  basehaz_line <- haz[atrisk$time_to_stop]

  cumhaz <- cumsum(haz)

  cumhaz_0_line <- cumhaz[atrisk$time_to_stop]
  cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]
  cumhaz_line <- (cumhaz_0_line - cumhaz_tstart)  * explp / newrisk

  chz_id_interval <- rowsum(cumhaz_line,
                            group = atrisk$id_interval,
                            reorder = TRUE) %>%
    as.data.frame()  %>%
    tibble::rownames_to_column() %>%
    tidyr::separate("rowname", into = c("id", "interval"), sep = "_", convert = TRUE) %>%
    arrange_("id", "interval")

  c_vecs <- split(chz_id_interval$V1, chz_id_interval$id)
#
#   browser()
#
#   poia<-nlm(f = dynfrail_fit, p =log(c(distribution$theta, distribution$lambda)),
#       dist = distribution$dist,
#       pvfm = distribution$pvfm, Y = Y, Xmat = X,
#       atrisk = atrisk, basehaz_line = basehaz_line,
#       mcox =list(coefficients = g, loglik = mcox$loglik),
#       c_vecs = c_vecs,
#       inner_control = control$inner_control)
#
#   optim(par = log(c(distribution$theta, distribution$lambda)), fn = dynfrail_fit,
#         dist = distribution$dist,
#         pvfm = distribution$pvfm, Y = Y, Xmat = X,
#         atrisk = atrisk, basehaz_line = basehaz_line,
#         mcox =list(coefficients = g, loglik = mcox$loglik),
#         c_vecs = c_vecs,
#         inner_control = control$inner_control)

  outer_m <- do.call(nlm,
                     args = c(list(f = dynfrail_fit, p = log(c(distribution$theta, distribution$lambda)),
                                   dist = distribution$dist,
                                   pvfm = distribution$pvfm, Y = Y, Xmat = X,
                                   atrisk = atrisk, basehaz_line = basehaz_line,
                                   mcox =list(coefficients = g, loglik = mcox$loglik),
                                   c_vecs = c_vecs,
                                   inner_control = control$inner_control), control$nlm_control))
#
#   outer_m <- poia
  inner_m <- dynfrail_fit(logfrailtypar = c(outer_m$estimate[1], outer_m$estimate[2]),
         dist = distribution$dist,
         pvfm = distribution$pvfm, Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
         mcox =list(coefficients = g, loglik = mcox$loglik),
         c_vecs = c_vecs,
         inner_control = control$inner_control,
         return_loglik = FALSE)

  res <- list(coefficients = inner_m$coef,
              hazard = inner_m$haz,
              imat = inner_m$Imat,
              logtheta = outer_m$estimate[1],
              loglambda = outer_m$estimate[2],
              # var_logtheta = NA,
              # ci_logtheta = NA,
              # var_loglambda = NA,
              # ci_loglambda = NA,
              frail_id = data.frame(id = df_dynfrail$id_,
                                    interval = df_dynfrail$interval_,
                                    frail = inner_m$frail,
                                    tstart = Y[,1],
                                    tstop = Y[,2],
                                    status = Y[,3]),
              # residuals = NA,
              tev = inner_m$tev,
              loglik = c(mcox$loglik[length(mcox$loglik)], -outer_m$minimum),
              formula = formula,
              distribution = distribution,
              control = control)
              # nobs = NA,
              # fitted = NA,
              # mf = NA,
              # mm = NA)

  attr(res, "call") <-  Call
  attr(res, "class") <- "dynfrail"


  res

  }
