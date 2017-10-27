#' Inner maximization of the log-likelihood
#'
#' @param logfrailtypar A vector containing the natural logarithm of the two parameters (\code{theta} - for the distribution, \code{lambda} - for the autocorrelation)
#' @param dist Argument of \code{\link{dynfrail_dist}}
#' @param pvfm Argument of \code{\link{dynfrail_dist}}
#' @param Y A \code{Surv} object obtained by splitting the original data at all the time points where the frailty process takes new values
#' @param Xmat A model matrix obtained by splitting the original data at all the time points where the frailty process takes new values
#' @param atrisk A list of various calculations that are used in the maximization process.
#' @param basehaz_line A vector with the baseline hazard estimate at each right hand side time point from \code{Y} (can be 0 for the others)
#' @param mcox An initial Cox model
#' @param c_vecs A list of the length equal to the number of clusters; each element contains a vector of the length of different values that \eqn{Z(t)} takes in that cluster.
#' Each element of this vector contains the sum of the cumulative hazards associated with that value of the frailty.
#' @param inner_control Argument of \code{\link{dynfrail_control}}
#' @param return_loglik Logical. If \code{TRUE}, then this just returns the log-likelihood, otherwise it returns also the estimates and information matrix
#'
#' @return The log-likelihood if \code{return_loglik = TRUE} or a list with the log-likelihood and estimates if \code{return_loglik = FALSE}.
#' @export
#'
#' @details This is an internal function that is used by \code{\link{dynfrail}}. It is not recommended to use this directly unless you know exactly what you are doing.
#' On the other hand, this might be useful if someone wants, for example, to use different maximizers, or to calculate the log-likelihood
#' at specific values of \code{theta, lambda}. Most of the input can be produced by \code{\link{dynfrail_prep}}.
#'
#' @examples
#' arglist1 <- dynfrail_prep(Surv(time, status) ~ rx + sex + cluster(litter),
#' data = rats)
#'
#' # using list() inside is because of the way that R converts lists and vectors
#' mod1 <- do.call(dynfrail_fit, c(logfrailtypar = list(log(c(0.5, 0.1))), arglist1))
dynfrail_fit <- function(logfrailtypar, #
                   dist,
                   pvfm,
                   Y,
                   Xmat,
                   atrisk,
                   basehaz_line,
                   mcox = list(),
                   c_vecs,
                   inner_control, # a list of some parameters
                   return_loglik = TRUE) {


  pars <- dist_to_pars(dist, logfrailtypar, pvfm)
  # browser()0.1
  if(length(Xmat)==0) {
    g_x <- matrix(rep(0, nrow(Y)),ncol = 1)
  } else {
    g_x <- t(mcox$coefficients %*% t(Xmat))
  }

  loglik_old = -Inf
  ncycles <- 0

  convergence <- FALSE
  while(!isTRUE(convergence)) {

    Estep <- lapply(seq_along(c_vecs), function(id) {
      Estep_id(events = atrisk$events_incluster[[id]], cvec = c_vecs[[id]],
               aalpha = pars$aalpha,
               ggamma = pars$ggamma, dist = pars$dist,
               pvfm = -1/2, times = atrisk$times_incluster[[id]], llambda = pars$llambda)
    })

    # log-likelihood
    llik_contrib <- sum(do.call(c, lapply(Estep, function(x) {
      log(abs(x[length(x) - 1])) + x[length(x)]
    })))

    loglik <- sum((log(basehaz_line) + g_x)[Y[,3] == 1]) + llik_contrib +
      sum(Y[,3]) - sum((atrisk$nevent * log(atrisk$nevent))[atrisk$nevent > 0])


    if(loglik < loglik_old - inner_control$lik_tol)
      stop(paste0("likelihood decrease of ", loglik - loglik_old ))

    if(abs(loglik - loglik_old) < inner_control$eps) break

    # print(loglik)
    #print(paste("beta", mcox$coefficients))
    loglik_old <- loglik


    # match the logz to the rows of the data frame
    # logz <- log(do.call(c, mapply(function(a, b) a[b],
    #                               lapply(Estep, function(x)
    #                                 -x[1:(length(x) - 2)] / x[length(x) - 1]),
    #                               atrisk$interval_incluster,
    #                               SIMPLIFY = FALSE
    # )))

    logz <- do.call(c,
    mapply(function(a,b) a[b],
           lapply(Estep, function(x) log(abs(x[1:(length(x) - 2)])) -
                                        log(abs(x[length(x) - 1]))),
           atrisk$interval_incluster,
           SIMPLIFY = FALSE)
    )

    mcox <- survival::agreg.fit(x = Xmat, y = Y, strata = NULL, offset = logz, init = NULL,
                                control = survival::coxph.control(), weights = NULL,
                                method = "breslow", rownames = NULL)

    # print(paste("logz[1:3]", logz[1], logz[2], logz[3]))

    if(length(Xmat)==0) {
      lp <- mcox$linear.predictors
      g_x <- t(matrix(rep(0, length(mcox$linear.predictors)), nrow = 1))
    } else {
      lp <- mcox$linear.predictors + as.numeric(t(mcox$coefficients) %*% mcox$means)
      g_x <- t(mcox$coefficients %*% t(Xmat))
    }


    explp <- exp(lp)
    #explp <- exp(mcox$linear.predictors)

    # newrisk <- exp(c(atrisk$x2 %*% mcox$coefficients) + 0)

    # Idea: nrisk has the sum of elp who leave later at every tstop
    # esum has the sum of elp who enter at every tstart
    # indx groups which esum is right after each nrisk;
    # the difference between the two is the sum of elp really at risk at that time point.


    nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))
    esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))


    nrisk <- nrisk - c(esum, 0,0)[atrisk$indx]
    haz <- atrisk$nevent/nrisk #  * newrisk
    # print(paste(haz[1], haz[2], haz[3], haz[4], haz[5]))


    cumhaz <- cumsum(haz)

    # baseline hazard for each tstop
    basehaz_line <- haz[atrisk$time_to_stop]
    cumhaz_0_line <- cumhaz[atrisk$time_to_stop]

    cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]
    cumhaz_line <- (cumhaz_0_line - cumhaz_tstart)  #* explp #  / newrisk

    chz_id_interval <- rowsum(cumhaz_line * exp(g_x),
                              group = atrisk$id_interval,
                              reorder = TRUE) %>%
      as.data.frame()  %>%
      tibble::rownames_to_column() %>%
      tidyr::separate("rowname", into = c("id", "interval"), sep = "_", convert = TRUE) %>%
      dplyr::arrange_("id", "interval")

    c_vecs <- split(chz_id_interval$V1, chz_id_interval$id)

    ncycles <- ncycles + 1
    if(ncycles > inner_control$maxit) {
      warning(paste("did not converge in ", inner_control$maxit," iterations." ))
      break
    }


    }

  # browser()

  if(isTRUE(return_loglik)) {
    # browser()
    if(isTRUE(inner_control$verbose))
      print(paste0("gamma: ",round(pars$ggamma, digits = 2),
                 " lambda: ",round(pars$llambda, digits = 2),
                 " loglik: ",round(loglik, digits = 4)))

    return(-loglik)
  }  # for when maximizing


  tev <- atrisk$time[haz > 0]
  haz_tev = haz[haz > 0]

  # if no SE, then return here

  nev_tp <- atrisk$nevent[atrisk$nevent!=0]

  z_elp <- exp(lp)
  elp = exp(lp)  / exp(logz)

  # browser()


  if(length(Xmat)>0) {
    x <- lapply(apply(Xmat, 1, list), function(x) x[[1]])
    x_z_elp <- Map(function(a,b) a*b, x, z_elp)
    x_z_elp_H0 <- Map(function(a,b,c) a*b*c, x, z_elp, cumhaz_line)
    x_elp_H0 <- Map(function(a,b,c) a*b*c, x, z_elp / exp(logz), cumhaz_line)

    xx <- lapply(x, function(x) x %*% t(x) )
    xx_z_elp_H0 <- Map(function(a,b, c) a * b * c, xx, z_elp, cumhaz_line)
    m_d2l_dgdg <- Reduce("+", xx_z_elp_H0)

    m_d2l_dhdg <-
      do.call(rbind,
              lapply(lapply(
                lapply(tev, function(tk) which(Y[,1] < tk & tk <= Y[,2])),
                function(x) x_z_elp[x]),
                function(...) Reduce("+", ...))
      )

  } else {
    m_d2l_dgdg <- NULL
    m_d2l_dhdg <- NULL
  }

  m_d2l_dhdh <- diag(nev_tp/haz_tev^2)


  Imat <- matrix(0, ncol(Xmat) + length(tev), ncol(Xmat) + length(tev))


  # Here the idea is to use indices instead of event times.
  # if tau1 = 0 and tau2 = 10, this means that the event time points for which that at risk
  # period stands for are 1, 2, 3, ... 10.
  # in C++ terms, this means positions 0 (from tau1) to 9 (so with tau2 we will always use < instead of <=)
  tau1 <- findInterval(Y[,1], tev)
  tau2 <- findInterval(Y[,2], tev, left.open = FALSE, rightmost.closed = FALSE)
  tau <- seq_along(tev)

  # todo: add cluster into the atrisk (or not who cares)
  cluster_id <- rep(1:length(atrisk$times_incluster), sapply(atrisk$interval_incluster, length))

  rows_tau <- lapply(split(data.frame(tau1, tau2), cluster_id), as.matrix)
  rows_elp <- split(elp, cluster_id)

  rows_x_elp_H0 <- lapply(split(as.data.frame(Xmat * elp * cumhaz_line), cluster_id), as.matrix)


  # within each individual, for each interval which lines are contained within that interval
  interval_rows <- lapply(atrisk$interval_incluster, function(x) {
    lapply(unique(x), function(y) which(y==x))
  })

  # Now to make sum calculations

  ez <- lapply(Estep, function(x)
    -x[1:(length(x) - 2)] / x[length(x) - 1])

  # browser()
  Iloss <- Vcov_adj(events_l = atrisk$events_incluster,
                       cvec_l = c_vecs,
                       aalpha = pars$aalpha,
                       ggamma = pars$ggamma, dist = pars$dist,
                       pvfm = -1/2, times_l = atrisk$times_incluster, llambda = pars$llambda,
                       elp_l = rows_elp,
                       xelph_l = rows_x_elp_H0,
                       tau_l = rows_tau,
                       interval_rows_l = interval_rows,
                       ez_l = ez,
                       n_times = length(tev),
                       n_covs = ncol(Xmat))

  # events_incluster is a vector with where the events are in the cluster, for the intervals of the frailty;
  # e.g. 1 1 3 means that there are two events in the FIRST frailty interval, one in 3.
  # this is relative to the intervals that exist in the cluster, not in all the intervals ever.
  # then cvec is a vector (c1, c2, c3) for example means that for 3 frailty intervals that is the summed up cumhaz
  # times_incluster assigns time points to the intervals (c1, c2, c3)
  # tau is sort of the Y but it delimits the rank of each time point instead of the time point itself
  # intervals_rows is a LIST: for each frailty interval it says which lines belong to that interval


  # id <- 163
  # Iloss_id <- Vcov_adj(events_l = atrisk$events_incluster[id],
  #                   cvec_l = c_vecs[id],
  #                   aalpha = pars$aalpha,
  #                   ggamma = pars$ggamma, dist = pars$dist,
  #                   pvfm = -1/2, times_l = atrisk$times_incluster[id], llambda = pars$llambda,
  #                   elp_l = rows_elp[id],
  #                   xelph_l = rows_x_elp_H0[id],
  #                   tau_l = rows_tau[id],
  #                   interval_rows_l = interval_rows[id],
  #                   ez_l = ez[id],
  #                   n_times = length(tev),
  #                   n_covs = ncol(Xmat))


if(!is.null(mcox$coefficients)) {
    Imat[1:length(mcox$coefficients), 1:length(mcox$coefficients)] <- m_d2l_dgdg - Iloss$betabeta
    Imat[1:length(mcox$coefficients), (length(mcox$coefficients)+1):nrow(Imat) ] <- t(m_d2l_dhdg) - t(Iloss$betalambda)
    Imat[(length(mcox$coefficients)+1):nrow(Imat), 1:length(mcox$coefficients) ] <- m_d2l_dhdg - Iloss$betalambda
  }

  # make it into matrix
  # https://stackoverflow.com/questions/37615790/indexing-the-elements-of-a-matrix-in-r

  Triangle1 <- function(k,n) {
    y <- -n
    r <- rep(0.0,n)
    t(vapply(1:n, function(x) {y <<- y+n+2L-x; c(rep(0L,x-1L),k[y:(y+n-x)])}, r))
  }

  cor_dh <- Triangle1(Iloss$lambdalambda, length(nev_tp))
  cor_dh[lower.tri(cor_dh)] <- t(cor_dh)[lower.tri(cor_dh)]

  Imat[(length(mcox$coefficients)+1):nrow(Imat), (length(mcox$coefficients)+1):nrow(Imat)] <- m_d2l_dhdh - cor_dh

  # browser()
  # solve(Imat) %>% diag %>% sqrt
  # try(solve(Imat))


  if(!isTRUE(return_loglik)) {
    return(list(loglik = loglik,
                tev = tev,
                haz = haz_tev,
                frail = exp(logz),
                coef = mcox$coefficients,
                Imat = Imat))
  }

}


