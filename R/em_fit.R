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
  # browser()
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
      warning(paste0("likelihood decrease of ", loglik - loglik_old ))

    if(abs(loglik - loglik_old) < inner_control$eps) break

    print(loglik)
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

    explp <- exp(mcox$linear.predictors)

    newrisk <- exp(c(atrisk$x2 %*% mcox$coefficients) + 0)

    # Idea: nrisk has the sum of elp who leave later at every tstop
    # esum has the sum of elp who enter at every tstart
    # indx groups which esum is right after each nrisk;
    # the difference between the two is the sum of elp really at risk at that time point.


    nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))
    esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))


    nrisk <- nrisk - c(esum, 0,0)[atrisk$indx]
    haz <- atrisk$nevent/nrisk * newrisk
    cumhaz <- cumsum(haz)

    # baseline hazard for each tstop
    basehaz_line <- haz[atrisk$time_to_stop]
    cumhaz_0_line <- cumhaz[atrisk$time_to_stop]

    cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]
    cumhaz_line <- (cumhaz_0_line - cumhaz_tstart)  * explp / newrisk

    chz_id_interval <- rowsum(cumhaz_line,
                              group = atrisk$id_interval,
                              reorder = TRUE) %>%
      as.data.frame()  %>%
      tibble::rownames_to_column() %>%
      tbl_df() %>%
      separate(rowname, into = c("id", "interval"), sep = "_", convert = TRUE) %>%
      arrange(id, interval)

    c_vecs <- split(chz_id_interval$V1, chz_id_interval$id)

    ncycles <- ncycles + 1
    if(ncycles > inner_control$maxit) {
      warning(paste("did not converge in ", inner_control$maxit," iterations." ))
      break
    }


    }


  if(isTRUE(return_loglik)) {
    if(isTRUE(inner_control$verbose)) print(paste("loglik = ",loglik))
    return(-loglik)
  }  # for when maximizing


}


