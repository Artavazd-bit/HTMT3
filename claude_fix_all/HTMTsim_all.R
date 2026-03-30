################################################################################
## simulation script for all sample sizes: 50, 200, 800, 3200
## including libpaths for running on HPC
## setup.R includes all functions used
## #############################################################################
.libPaths(c("/home/jab49wd/R-projects/R/library", .libPaths()))

library(lavaan)
library(doParallel)
library(foreach)
library(dplyr)

# Requires setup.R in the working directory; deploy both files to the same HPC directory
source("setup.R")

# Log session info
writeLines(capture.output(sessionInfo()), "sessioninfo_HTMTsim_all.txt")

##################### Simulation Design ########################################
# simModels is defined in setup.R; switch to simModels_tau for tau-equivalent models
nkernel <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "32"))
nobs <- c(50, 200, 800, 3200, 12800)
simrunstotal <- 1000
bootruns <- 1000
alphavec <- c(0.01, 0.05, 0.1)
# normal: no skewness, no kurtosis, nonnormal: specified skewness and kurtosis
# values
disttable <- data.frame(name = c("normal", "nonnormal"))
disttable$skewness <- list(NULL, c(0.7, 0.9, 1.2, 1.5, 0.8, 1.3))
disttable$kurtosis <- list(NULL, c(3.5, 3.6, 3.7, 3.5, 3.6, 3.7))
######################### Monte Carlo Simulation ###############################
cl <- parallel::makeCluster(nkernel, outfile = "errorcluster.txt")
on.exit(parallel::stopCluster(cl), add = TRUE)
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/jab49wd/R-projects/R/library"))

# iteration order: simModels -> samplesize -> 10000 simruns ->
# datatype(normal, nonnormal)
# for each simrun a deterministic seed is used that encodes the condition
# (model, sample size, replication) to ensure independence across conditions.
# the wrapper(...) function creates the percentile, asymptotic and BCa
# confidence intervals for multiple confidence levels.

# Incremental saving with explicit outer loops
# Clean results directory to prevent stale .rds files from previous runs
dir.create("results", showWarnings = FALSE)
old_rds <- list.files("results", pattern = "partial_.*\\.rds", full.names = TRUE)
if (length(old_rds) > 0) file.remove(old_rds)

for (jj in 1:nrow(simModels)) {
  for (n in nobs) {

    # Set up variables that depend on jj and n for export to workers
    current_model <- simModels$model[jj]
    current_correlation <- simModels$correlation[jj]

    partial <- foreach(sim_runs = 1:simrunstotal, .combine = "rbind",
                       .packages = c("lavaan", "foreach", "dplyr"),
                       .export = c("wrapper", "deltamethod", "derivhtmt", "calchtmt",
                                   "calcovcov", "calcovcor", "bootbca", "bootstrap",
                                   "jacknife", "simModels", "model_est", "bootruns",
                                   "alphavec", "disttable", "simrunstotal",
                                   "current_model", "current_correlation")) %dopar%
    {
      tryCatch({
        # Deterministic seed encoding model (jj), sample size (n), and replication
        seed <- sim_runs + (jj - 1) * 10000 + n * 100000
        # Set RNG state so bootstrap/jackknife resampling is also reproducible
        set.seed(seed)
        correlation <- current_correlation
        temp <- foreach(distn = 1:nrow(disttable), .packages = c("lavaan", "foreach", "dplyr"), .combine = "rbind") %do%
          {
            data <- lavaan::simulateData(model = current_model,
                                         sample.nobs = n,
                                         skewness = disttable$skewness[[distn]],
                                         kurtosis = disttable$kurtosis[[distn]],
                                         seed = seed,
                                         empirical = FALSE,
                                         return.type = "data.frame"
            )

            simuresults <- wrapper(data = data, model = model_est, latent1 = "xi_1", latent2 = "xi_2", alpha = alphavec, scale = FALSE, htmt2 = FALSE, nboot = bootruns)
            res <- foreach(type = c("delta", "boot", "bcaboot"), .combine = "rbind") %:%
              foreach(indexgamma = 1:length(alphavec), .combine = "rbind") %do%
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
                                  , sim_runs
                                  , datatype = disttable$name[distn]
                                  , method = type
                                  , alpha = alphavec[indexgamma]
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
        temp
      }, error = function(e) {
        cat(paste("Error in sim_runs", sim_runs, ":", conditionMessage(e), "\n"))
        data.frame(
          correlation = NA, n = NA, HTMT = NA, seHTMT = NA,
          sim_runs = sim_runs, datatype = NA, method = NA, alpha = NA,
          lowerbound = NA, upperbound = NA,
          coveragecorr = NA, coverageone = NA,
          time = NA, seed = NA, missing = NA
        )
      })
    }

    saveRDS(partial, paste0("results/partial_model", jj, "_n", n, ".rds"))
    cat(paste("Completed model", jj, "n =", n, "\n"))
  }
}

# Assemble final CSV from saved partials
ergebinger <- do.call(rbind, lapply(
  list.files("results", pattern = "partial_.*\\.rds", full.names = TRUE), readRDS
))

# stopCluster also handled by on.exit() as safety net
parallel::stopCluster(cl)
write.csv2(ergebinger, "HTMTsimresults_all.csv")
