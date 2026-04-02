library(foreach)
library(lavaan)
library(dplyr)
source("arraysimulation/setup.R")

conditions <- expand.grid(rep_batch = 1:10,
                          model = 1:nrow(simModels), 
                          n = c(25, 50, 100, 200, 400, 800, 1600, 3200, 6400),
                          datatype = c(1,2)
)

saveRDS(conditions, "arraysimulation/conditions.rds")
