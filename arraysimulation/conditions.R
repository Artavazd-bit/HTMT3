
source("arraysimulation/setup.R")

conditions <- expand.grid(rep_batch = 1:10,
                          model = 1:nrow(simModels), 
                          n = c(50, 100, 250, 1000, 6000),
                          stringsAsFactors = FALSE,
                          datatype = c(1,2)
)

saveRDS(conditions, "arraysimulation/conditions.rds")
