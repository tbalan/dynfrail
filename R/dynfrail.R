dynfrail <- function(.data,
                     .formula,
                     .distribution = dynfrail_distribution(),
                     .control = dynfrail_control()) {


  Call <- match.call()


  if(missing(.formula)  | missing(.data)) stop("Missing arguments")

  cluster <- function(x) x

  mf_tmp <- model.frame(.formula, .data)


  tev_unique_ord <- sort(unique(mf_tmp[[1]][,dim(mf_tmp[[1]])[2] - 1][mf_tmp[[1]][,dim(mf_tmp[[1]])[2]] == 1]))

  if(!is.null(.distribution$n_ints)) {
    cut <- findInterval(x = quantile(tev_unique_ord, probs = seq(from = 0, to = 1, length.out = n_ints + 2))[-c(1, n_ints + 1)],
      vec = tev_unique_ord)

  } else
    if(!is.null(.distribution$times)) cut <- .distribution$times else
      cut <- tev_unique_ord
      # then cut is the unique event time points


  .df_dynfrail <- survSplit(.formula, data = .data, cut = cut)
  #name1 <- names(.df_dynfrail[grep("cluster", names(.df_dynfrail))])

  .df_dynfrail$id <- .df_dynfrail$`cluster(id)`
  # this has to be renamed. figure out how!
  mf <- model.frame(.formula, .df_dynfrail)

  # Identify the cluster and the ID column // this might have to happen at some point earlier
  pos_cluster <- grep("cluster", names(mf))
  if(length(pos_cluster) != 1) stop("misspecified or non-specified cluster")
  id <- mf[[pos_cluster]]

  Y <- mf[[1]]

  # check Y really. in this case it can't be with 2 columns and it should not be!
  if(!inherits(Y, "Surv")) stop("left hand side not a survival object")
  if(ncol(Y) != 3) {
    Y <- Surv(rep(0, nrow(Y)), Y[,1], Y[,2])
  }
  if(attr(Y, "type") != "counting") stop("use Surv(tstart, tstop, status)")



  # get the model matrix
  X1 <- model.matrix(.formula, .df_dynfrail)
  # this is necessary because when factors have more levels, pos_cluster doesn't correspond any more
  pos_cluster_X1 <- grep("cluster", colnames(X1))
  pos_terminal_X1 <- grep("terminal", colnames(X1))
  X <- X1[,-c(1, pos_cluster_X1, pos_terminal_X1), drop=FALSE]
  # note: X has no attributes, in coxph it does.

  # some stuff for creating the C vector, is used all along.
  # mcox also works with empty matrices, but also with NULL as x.
  mcox <- survival::agreg.fit(x = X, y = Y, strata = NULL, offset = NULL, init = NULL,
                              control = survival::coxph.control(),
                              weights = NULL, method = "breslow", rownames = NULL)


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

  explp <- exp(mcox$linear.predictors) # these are with centered covariates

  nev_id <- rowsum(Y[,3], id) # nevent per id or am I going crazy



  nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))
  esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))

  # the stuff that won't change
  death <- (Y[, ncol(Y)] == 1)
  nevent <- as.vector(rowsum(1 * death, Y[, ncol(Y) - 1])) # per time point
  time <- sort(unique(Y[,2])) # unique tstops

  # this gives the next entry time for each unique tstop (not only event)
  etime <- c(0, sort(unique(Y[, 1])),  max(Y[, 1]) + min(diff(time)))
  indx <- findInterval(time, etime, left.open = TRUE) # left.open  = TRUE is very important

  # this gives for every tstart (line variable) after which event time did it come
  # indx2 <- findInterval(Y[,1], time, left.open = FALSE, rightmost.closed = TRUE)
  indx2 <- findInterval(Y[,1], time)

  time_to_stop <- match(Y[,2], time)
  order_id <- findInterval(id, unique(id))

  atrisk <- list(death = death, nevent = nevent, nev_id = nev_id,
                 order_id = order_id, time = time, indx = indx, indx2 = indx2,
                 time_to_stop = time_to_stop)

  nrisk <- nrisk - c(esum, 0,0)[indx]

  haz <- nevent/nrisk * newrisk


  basehaz_line <- haz[atrisk$time_to_stop]

  cumhaz <- cumsum(haz)

  cumhaz_0_line <- cumhaz[atrisk$time_to_stop]
  cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]
  cumhaz_line <- (cumhaz_0_line - cumhaz_tstart)  * explp / newrisk

  Cvec <- rowsum(cumhaz_line, atrisk$order_id)


  # now here it goes: with these ingredients we are ready to go to the EM!

  # first: a one fit for fixed theta and lambda
}
