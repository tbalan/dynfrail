em_fit <- function(logfrailtypar, # a vector of two parameters (theta - for the distribution, llambda - for the distance)
                   dist,
                   pvfm,
                   Y,
                   Xmat,
                   atrisk,
                   basehaz_line,
                   mcox = list(),
                   c_vecs,
                   times,
                   inner_control, # a list of some parameters
                   se = FALSE, # whether standard errors should be calculated
                   return_loglik = TRUE) {


  pars <- dist_to_pars(dist, logfrailtypar, pvfm)
  # browser()0.1
  if(length(Xmat)==0) {
    g_x <- matrix(rep(0, nrow(Y)),ncol = 1)
  } else {
    g_x <- t(mcox$coefficients %*% t(Xmat))
  }

  # Here a check whether logfrailtypar is at the edge of the parameter space
  # if(logfrailtypar > inner_control$lower_tol) {
  #   #message("Frailty parameter very large, frailty variance close to 0")
  #   loglik <- mcox$loglik[length(mcox$loglik)]
  #   # loglik <- sum((log(basehaz_line) + g_x)[Y[,3] == 1]) +
  #   #    sum(Y[,3]) - sum(nev_tp * log(nev_tp))
  #
  #   if(isTRUE(return_loglik)) {
  #     if(isTRUE(inner_control$verbose)) print(paste("loglik = ",loglik))
  #     return(-loglik)
  #   }
  #
  # }


  loglik_old = -Inf
  ncycles <- 0
  # browser()

  convergence <- FALSE
  while(!isTRUE(convergence)) {

    Estep <- lapply(seq_along(c_vecs), function(id) {
      Estep_id(events = atrisk$events_incluster[[id]], cvec = c_vecs[[id]],
               aalpha = pars$aalpha,
               ggamma = pars$ggamma, dist = 0,
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
    logz <- log(do.call(c, mapply(function(a, b) a[b],
                                  lapply(Estep, function(x)
                                    -x[1:(length(x) - 2)] / x[length(x) - 1]),
                                  atrisk$interval_incluster,
                                  SIMPLIFY = FALSE
    )))


    mcox <- survival::agreg.fit(x = Xmat, y = Y, strata = NULL, offset = logz, init = NULL,
                                control = survival::coxph.control(), weights = NULL,
                                method = "breslow", rownames = NULL)


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
      tidyr::separate(rowname, into = c("id", "interval"), sep = "_", convert = TRUE) %>%
      dplyr::arrange(id, interval)

    c_vecs <- split(chz_id_interval$V1, chz_id_interval$id)

    ncycles <- ncycles + 1
    if(ncycles > inner_control$maxit) {
      warning(paste("did not converge in ", inner_control$maxit," iterations." ))
      break
    }


    }


  if(isTRUE(return_loglik)) {
    # browser()

    print(paste0("ggamma: ",round(pars$ggamma, digits = 2),
                 " llambda: ",round(pars$llambda, digits = 2),
                 " loglik: ",round(loglik, digits = 3)))

    if(isTRUE(inner_control$verbose)) print(paste("loglik = ",loglik))
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

  # if(!is.null(mcox$coefficients)) {
  #   Imat[1:length(mcox$coefficients), 1:length(mcox$coefficients)] <- m_d2l_dgdg
  #   Imat[1:length(mcox$coefficients), (length(mcox$coefficients)+1):nrow(Imat) ] <- t(m_d2l_dhdg)
  #   Imat[(length(mcox$coefficients)+1):nrow(Imat), 1:length(mcox$coefficients) ] <- m_d2l_dhdg
  # }
  #
  # Imat[(length(mcox$coefficients)+1):nrow(Imat), (length(mcox$coefficients)+1):nrow(Imat)] <- m_d2l_dhdh

  # This is d/dg

  # dl1_dg <- apply(Xmat * Y[,3], 2, sum)
  #
  # dl2_dg <- apply(Xmat  * z_elp * cumhaz_line, 2, sum)
  #
  # dl1_dl <- sum(nev_tp / haz_tev)

  # all this stuff is 0 in the end, isn't it?

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

  ez <- lapply(Estep, function(x) -x[1:(length(x) - 2)] / x[length(x) - 1])

  Iloss <- Vcov_adj(events = atrisk$events_incluster,
                       cvec = c_vecs,
                       aalpha = pars$aalpha,
                       ggamma = pars$ggamma, dist = 0,
                       pvfm = -1/2, times = atrisk$times_incluster, llambda = pars$llambda,
                       elp = rows_elp,
                       xelph = rows_x_elp_H0,
                       tau = rows_tau,
                       interval_rows = interval_rows,
                       ez = ez,
                       n_times = length(tev),
                       n_covs = ncol(Xmat))


  if(!is.null(mcox$coefficients)) {
    Imat[1:length(mcox$coefficients), 1:length(mcox$coefficients)] <- m_d2l_dgdg + Iloss$betabeta
    Imat[1:length(mcox$coefficients), (length(mcox$coefficients)+1):nrow(Imat) ] <- t(m_d2l_dhdg) + t(Iloss$betalambda)
    Imat[(length(mcox$coefficients)+1):nrow(Imat), 1:length(mcox$coefficients) ] <- m_d2l_dhdg + Iloss$betalambda
  }

  # make it into matrix
  # https://stackoverflow.com/questions/37615790/indexing-the-elements-of-a-matrix-in-r

  Triangle1 <- function(k,n) {
    y <- -n
    r <- rep(0.0,n)
    t(vapply(1:n, function(x) {y <<- y+n+2L-x; c(rep(0L,x-1L),k[y:(y+n-x)])}, r))
  }

  Imat[(length(mcox$coefficients)+1):nrow(Imat), (length(mcox$coefficients)+1):nrow(Imat)] <- m_d2l_dhdh + Triangle1(Iloss$lambdalambda, length(nev_tp))

  se <- try(sqrt(diag(solve(Imat))))


  if(!isTRUE(return_loglik)) {
    return(list(mcox = mcox, frail = exp(logz), cumhaz = cumhaz, se = se))

  }



}


