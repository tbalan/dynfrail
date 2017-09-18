data(rats)


head(cgd)
library(survival)
library(tidyverse)
source('~/dynfrail/R/dynfrail_arguments.R')

dynfrail(formula = Surv(tstart, tstop, status) ~ sex + treat + cluster(id),
         data = cgd
         )
