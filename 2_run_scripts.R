source("config.R")


print("age")
source("2_model_run_MTGP_age.R")
rm(list = ls())

print("hisp")
source("2_model_run_MTGP_hisp.R")
rm(list = ls())

print("race")
source("2_model_run_MTGP_race.R")
rm(list = ls())

print("sex")
source("2_model_run_MTGP_sex.R")
rm(list = ls())

print("total")
source("2_model_run_MTGP_total.R")
rm(list = ls())

print("type")
source("2_model_run_MTGP_type.R")
rm(list = ls())
