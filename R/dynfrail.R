dynfrail <- function(formula, data,
                     distribution = dynfrail_distribution(),
                     control = dynfrail_control())
{


  Call <- match.call()

  # Here checking inputs
  if(missing(formula)  | missing(data)) stop("Missing arguments")

  cluster <- function(x) x

  Y <- model.frame(formula, data)[[1]]

  if(ncol(Y) != 3)
    stop("Provide (tstart, tstop, status) format")

  tev_unique_ord <- sort(unique(Y[,2][Y[,3]==1]))

  if(!is.null(distribution$n_ints)) {
    cut <- findInterval(x = quantile(tev_unique_ord, probs = seq(from = 0, to = 1, length.out = n_ints + 2))[-c(1, n_ints + 1)],
      vec = tev_unique_ord)

  } else
    if(!is.null(distribution$times)) cut <- .distribution$times else
      cut <- tev_unique_ord

  df_dynfrail <- survSplit(formula, data = data, cut = cut, episode = "interval_")
  names(df_dynfrail)[grep("cluster", names(df_dynfrail))] <- "id_"



  terms2 <- drop.terms(terms(formula), drop = 3, keep.response = TRUE)

  mf <- model.frame(terms2, df_dynfrail)

  Y <- mf[[1]]
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

  nrisk <- nrisk - c(esum, 0,0)[indx]

  haz <- nevent/nrisk * newrisk

  basehaz_line <- haz[atrisk$time_to_stop]

  cumhaz <- cumsum(haz)

  cumhaz_0_line <- cumhaz[atrisk$time_to_stop]
  cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]
  cumhaz_line <- (cumhaz_0_line - cumhaz_tstart)  * explp / newrisk


  # Now to build up the c vectors
  id_interval <- paste0(df_dynfrail$id_, "_",df_dynfrail$interval_)

  c_vecs <- split(rowsum(cumhaz_line, group = id_interval), df_dynfrail$id_)
  delta <- split(rowsum(Y[,3], group = id_interval), df_dynfrail$id_)

  intervals <- split(df_dynfrail$interval_, df_dynfrail$id_)
  tmp <- lapply(delta, function(x) {
    rep(1:length(x), x)
  })

  # Now: the thing is that there are some missing but we don't care about that

    for(i in 1:length(c_vecs)) {
    c_vecs[[1]]

    SdivideC(tmp = rep(1:length(delta[[i]], delta[[i]] )))
  }

  # this has the chz for each tau interval from each person


  # next steps in C++:
  # take this and then loop over individuals
  # determine the vector of events that we need to
  # determine the maximum number of intervals for the taus (we have that probably

  # then select what's relevant for one individual
  # then calculate the expansion of the lambdas

  Cvec <- rowsum(cumhaz_line, atrisk$order_id)


  ## Input for E step:
  # -
  # now here it goes:

  # first: a one fit for fixed theta and lambda
}
