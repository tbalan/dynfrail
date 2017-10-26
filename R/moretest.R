# library(parfm)
# library(dplyr)
# data(asthma)
#
# small_asthma <-
#   asthma %>%
#   group_by(Patid) %>%
#   mutate(linenr = 1:n()) %>%
#   filter(linenr <= 3)
#
# library(dynfrail)
#
# # m1 <- dynfrail(Surv(Begin, End, Status) ~ Drug + cluster(Patid), data = small_asthma,
# #                distribution = dynfrail_dist(n_ints = 0, theta = 2))
# # m1
#
# m2_ig <- dynfrail(Surv(Begin, End, Status) ~ Drug + cluster(Patid), data = small_asthma,
#                   distribution = dynfrail_dist(dist = "pvf", n_ints = 2),
#                   control = dynfrail_control(inner_control = list(verbose = TRUE)))
# m2_ig
#
# m2_gam <- dynfrail(Surv(Begin, End, Status) ~ Drug + cluster(Patid), data = small_asthma,
#                   distribution = dynfrail_dist(dist = "gamma", n_ints = 2))
# m2_gam
#
#
# m3_gam <- dynfrail(Surv(Begin, End, Status) ~ Drug + cluster(Patid), data = small_asthma,
#                    distribution = dynfrail_dist(dist = "gamma", n_ints = 3), control = dynfrail_control(inner_control = list(verbose = TRUE)))
# m3_gam
#
#
# m3_ig <- dynfrail(Surv(Begin, End, Status) ~ Drug + cluster(Patid), data = small_asthma,
#                    distribution = dynfrail_dist(dist = "pvf", n_ints = 3))
# m3_ig
#
# m5_ig <- dynfrail(Surv(Begin, End, Status) ~ Drug + cluster(Patid), data = small_asthma,
#                   distribution = dynfrail_dist(dist = "pvf", n_ints = 4), control = dynfrail_control(inner_control = list(verbose = TRUE)))
#
#
#
# #
# # library(frailtyEM)
# # m1_sf <- emfrail(Surv(Begin, End, Status) ~ Drug + cluster(Patid), data = small_asthma)
# # summary(m1_sf)
# #
# # m1_cph <- coxph(Surv(Begin, End, Status) ~ Drug + frailty(Patid), data = small_asthma, ties = "breslow")
# # summary(m1_sf)
# #
# # m1_cph$coefficients
# # m1_sf$coefficients
# # m1$coefficients
# #
# #
# # a <- m1_cph$history$`frailty(Patid)`$c.logli
# # b <- m1_sf$loglik[2]
# # c <- m1$loglik[2]
# #
# # plot(c(a,b,c))
# #
# # b>c
