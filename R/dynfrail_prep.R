#' Preparation of the input for \code{dynfrail_fit}
#'
#' @param formula A formula that contains on the left hand side an object of the type \code{Surv}
#' and on the right hand side a \code{+cluster(id)} statement.
#' @param data A data frame in which the formula argument can be evaluated
#' @param distribution An object as created by \code{\link{dynfrail_dist}}
#' @param control An object as created by \code{\link{dynfrail_control}}
#' @param ... Other arguments, currently used to warn about deprecated argument names
#'
#' @return A list with what is needed to be used with \code{\link{dynfrail_fit}}
#' @export
#'
#' @details This is an internal function of \code{\link{dynfrail}} thath actually does before going to the inner maximization, except for the starting values.
#' The input is identical to that from \code{dynfrail}
#' A scenario where this would be useful would be to make these calculations and then things would be passed on to \code{\link{dynfrail_fit}}.
#'
#' @examples
#' arglist1 <- dynfrail_prep(Surv(time, status) ~ rx + sex + cluster(litter),
#' data = rats)
#'
#' @seealso \code{\link{dynfrail}}, \code{\link{dynfrail_fit}}
dynfrail_prep <- function(formula, data,
                          distribution = dynfrail_dist(),
                          control = dynfrail_control(),
                          ...) {


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
                                   length.out = distribution$n_ints + 2))


    cut <- tev_unique_ord[findInterval(x = quants[-c(1, length(quants))],
                                       vec = tev_unique_ord)]

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

  atrisk$events_incluster <- lapply(delta, function(x) {
    rep(1:length(x), x)
  })

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

  list(
       dist = distribution$dist,
       pvfm = distribution$pvfm,
       Y = Y, Xmat = X,
       atrisk = atrisk, basehaz_line = basehaz_line,
       mcox =list(coefficients = g, loglik = mcox$loglik),
       c_vecs = c_vecs,
       inner_control = control$inner_control
       )
}
