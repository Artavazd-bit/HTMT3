library(lavaan)
library(covsim)
library(rvinecopulib)
library(foreach)

source("arraysimulation/setup.R")

marginsxi <- lapply(X = c(4, 5, 6, 4, 5, 6),
                    function(X) list(distr = "chisq", df = X))

for (model_id in 1:nrow(simModels)) {
  current_model <- simModels$model[model_id]
  fit <- sem(current_model, do.fit = FALSE)
  popcov <- lavInspect(fit, "implied")$cov

  set.seed(model_id)
  vine <- vita(marginsxi,
               sigma.target = popcov,
               Nmax = 1000000)

  outfile <- paste0("arraysimulation/vines/vine_model_", model_id, ".rds")
  saveRDS(vine, outfile)
  cat("Saved vine for model", model_id, "\n")
}
