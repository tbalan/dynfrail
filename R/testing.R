data(rats)


head(cgd)
library(survival)
library(tidyverse)
source('~/dynfrail/R/dynfrail_arguments.R')
source("R/dynfrail.R")
Rcpp::sourceCpp("estep_new.cpp")

dynfrail(formula = Surv(tstart, tstop, status) ~ sex + treat + cluster(id),
         data = cgd
         )

data(rats)
head(rats)
dynfrail(formula = Surv(rep(0, nrow(rats)), time, status) ~ sex + rx + cluster(litter),
         data = rats
)
sums(1:5, 2, 4)
