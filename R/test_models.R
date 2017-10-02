# library(tidyverse)
# library(dynfrail)
# data(asthma)
# small_asthma <- asthma %>%
#   group_by(Patid) %>%
#   mutate(rn = row_number()) %>% ungroup() %>%
#   filter(rn <= 4) %>% mutate(Begin = Begin / 10, End = End / 10) %>%
#   mutate(y=rnorm(n()))
#
#
# library(frailtyEM)
# emfrail(Surv(Begin, End, Status) ~ cluster(Patid) + Drug, data = small_asthma,
#         distribution = emfrail_dist(theta = 2),
#         control = emfrail_control(opt_fit = FALSE))
#
# #
# # log-likelihood: -4023.475
# # theta: 2.461172
#
#
# mod1 <- dynfrail(Surv(Begin, End, Status) ~ cluster(Patid) + Drug + y, data = small_asthma,
#          distribution = dynfrail_distribution(n_ints = 3, theta = 2),
#          control = dynfrail_control(inner_control = list(maxit = 10)))
#
# mod1 <- dynfrail(Surv(Begin, End, Status) ~ cluster(Patid) + Drug, data = small_asthma,
#                  distribution = dynfrail_distribution(dist = "stable", n_ints = 3, theta = 2),
#                  control = dynfrail_control(inner_control = list(maxit = 10)))
#
# mod1 <- dynfrail(Surv(Begin, End, Status) ~ cluster(Patid) + Drug, data = small_asthma,
#                  distribution = dynfrail_distribution(dist = "pvf", n_ints = 3, theta = 2),
#                  control = dynfrail_control(inner_control = list(maxit = 10)))
#
# mod1$outer_m$minimum
#
#
