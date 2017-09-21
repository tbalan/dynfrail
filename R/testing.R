
library(survival)
library(tidyverse)
source('~/dynfrail/R/dynfrail_arguments.R')
source("R/dynfrail_aux.R")
source("R/dynfrail.R")
source('~/dynfrail/R/em_fit.R')
Rcpp::sourceCpp("estep_new.cpp")


asthma <- read.table("asthma.txt")
head(asthma)

small_asthma <- asthma %>%
  group_by(Patid) %>%
  mutate(rn = row_number()) %>% ungroup() %>%
  filter(rn <= 4) %>% mutate(Begin = Begin / 10, End = End / 10)

nrow(small_asthma)

dynfrail(Surv(Begin, End, Status) ~ Drug + cluster(Patid),
         data = small_asthma,
         distribution = dynfrail_distribution(times =
                                                c(6.01, 12.02, 18.03, 24.04, 30.05, 36.06, 42.07, 48.08, 54.09)))



dynfrail(formula = Surv(tstart, tstop, status) ~ sex + treat + cluster(id),
         data = cgd, distribution = dynfrail_distribution(n_ints = 3)
         )

data(rats)
head(rats)
rats <- rats %>%
  mutate(tstart = 0)
dynfrail(formula = Surv(tstart, time, status) ~ sex + rx + cluster(litter),
         data = rats,distribution = dynfrail_distribution(theta = 2, lambda = 0.1)
)

library(frailtyEM)
emfrail(formula = Surv(tstart, time, status) ~ sex + rx + cluster(litter),
         data = rats
)


sums(1:5, 2, 4)
