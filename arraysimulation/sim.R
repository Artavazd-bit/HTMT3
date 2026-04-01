.libPaths(c("/home/jab49wd/R-projects/R/library", .libPaths()))

library(lavaan)
library(dplyr)
library(foreach)
library(covsim)
library(rvinecopulib)

source("code/setup.R")
conditions <- readRDS("code/conditions.rds")
n_batches <- max(conditions$rep_batch)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cond <- conditions[task_id,]

total_reps <- 1000
number_bootruns <- 1000
alpha_vec <- c(0.05, 0.1, 0.01)

rep_batch <- cond$rep_batch
reps_per_batch <- total_reps / n_batches
start_rep <- (cond$rep_batch - 1) * reps_per_batch + 1  
end_rep <- (cond$rep_batch) * reps_per_batch


n <- cond$n
n_id <- match(n, c(50, 100, 250, 1000, 6000))
model_id <- cond$model
data_id <- cond$datatype
current_model <- simModels$model[model_id]
correlation <- simModels$correlation[model_id]
latent1 <- "xi_1"
latent2 <- "xi_2"


if(data_id == 1){
  datatype <- "normal"
}else if(data_id == 2){
  datatype <- "nonnormal"
}

if(data_id == 2){
  fit <- sem(current_model, do.fit = FALSE)
  popcov <- lavInspect(fit, "implied")$cov
  marginsxi <- lapply(X = c(4, 5, 6, 4, 5, 6), 
                      function(X) list(distr = "chisq", df = X))
  seed_vine <- as.integer(model_id * 100 + n_id)
  set.seed(seed_vine)
  vine <- vita(marginsxi, 
               sigma.target = popcov, 
               Nmax = 1000000)
  model_df <- lavaanify(current_model)
  listind1 <- list(model_df$rhs[model_df$lhs == latent1 & model_df$op == "=~"])
  listind2 <- list(model_df$rhs[model_df$lhs == latent2 & model_df$op == "=~"])
  all_indicators <- unlist(list(listind1, listind2)) 
  
  columnames <- all_indicators
}


partial <- foreach(r = start_rep : end_rep, .combine = "rbind")%do%
      {
        seed <- r + (model_id) * 10000 + (data_id) * 100000 + (n_id) * 1000000
        seed <- as.integer(seed)
        if(data_id == 1){
        data <- lavaan::simulateData(model = current_model,
                               sample.nobs = n,
                               seed = seed,
                               skewness = NULL, 
                               kurtosis = NULL, 
                               empirical = FALSE,
                               return.type = "data.frame"
        )
        }else if(data_id == 2){
          set.seed(seed)
          data <- rvinecopulib::rvine(n = n,
                                      vine = vine)
          data <- as.data.frame(data)
          colnames(data) <- columnames
        }
        simuresults <- wrapper(data = data,
                        model = model_est,
                        latent1 = latent1,
                        latent2 = latent2,
                        alpha = alpha_vec,
                        scale = FALSE,
                        htmt2 = FALSE,
                        nboot = number_bootruns)
        
        res <- foreach(type = c("delta", "boot", "bcaboot", "bcboot"), .combine = "rbind")%:%
                foreach(indexgamma = 1:length(alpha_vec), .combine = "rbind") %do%
                {
                  HTMT <- simuresults$delta$HTMT
                  seHTMT <- simuresults[[type]]$se
                  lowerbound <- unname(simuresults[[type]]$lowerbound[indexgamma])
                  upperbound <- unname(simuresults[[type]]$upperbound[indexgamma])
                  missing <- unname(simuresults[[type]]$missing)
                  time <- simuresults[[type]]$time
                  row <- data.frame(correlation = correlation
                                    , n = n
                                    , HTMT = HTMT
                                    , seHTMT = seHTMT
                                    , simruns = r
                                    , batch = rep_batch
                                    , datatype = datatype
                                    , method = type
                                    , alpha = alpha_vec[indexgamma]
                                    , lowerbound = lowerbound
                                    , upperbound =  upperbound
                                    , coveragecorr = lowerbound <= correlation & correlation <=  upperbound
                                    , coverageone =  lowerbound <= 1 & 1 <= upperbound
                                    , time = time
                                    , seed = seed
                                    , missing = missing
                  )
                  row
                }
          res
      }

outfile <- paste0("~/simulationhtmt/results/model_", cond$model, "_sample_size_",
                  cond$n, "_batch_", cond$rep_batch, "_datatype_", data_id, ".rds")
tmpfile <- paste0(outfile, ".tmp")
saveRDS(partial, tmpfile, version = 2)
file.rename(tmpfile, outfile)
        
        
        
