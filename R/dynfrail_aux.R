#' dist_to_pars
#'
#' @param dist One of gamma, stable, pvf
#' @param logfrailtypar log theta and log llambda
#' @param pvfm The pvfm
#'
#' @keywords internal
#' @return A list with 4 elements: aalpha, ggamma (the parameters of the Laplace transform), llambda (distance parameter) and dist_id.
dist_to_pars <- function(dist, logfrailtypar, pvfm) {

  if (dist == "gamma") {
    aalpha <- ggamma <- exp(logfrailtypar[1])
    dist_id <- 0L
  }

  # if (dist == "stable") {
  #     theta <- exp(logfrailtypar) + 1  # so theta >1
  #     bbeta <- 1 - 1/theta  # so bbeta in (0,1), that's what's important
  #     alpha <- theta / (theta - 1)  # alpha = 1/beta for scaling
  #     dist_id <- 1L
  # }

  if (dist == "stable") {
    # theta <- exp(logfrailtypar) + 1 # so theta >1
    # bbeta <- 1 - 1/theta
    aalpha <- 1
    #bbeta <- 1 - exp(logfrailtypar) / (exp(logfrailtypar) + 1)
    ggamma <- exp(logfrailtypar[1]) / (exp(logfrailtypar[1]) + 1)
    dist_id <- 1L
  }

  if (dist == "pvf") {
    aalpha <- abs((pvfm + 1)/pvfm * exp(logfrailtypar[1]))
    ggamma <- (pvfm + 1) * exp(logfrailtypar[1])
    dist_id <- 2L
  }

  list(aalpha = aalpha, ggamma = ggamma, llambda = exp(logfrailtypar[2]), dist = dist_id)
}
