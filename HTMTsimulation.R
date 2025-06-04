.libPaths(c("/home/jab49wd/R-projects/R/library", .libPaths()))

library(lavaan)
library(doParallel)
library(foreach)
library(dplyr)

source("setup.R")

##################### Simulation Design ########################################
simModels <- simModels
nkernel <- 32
nobs <- c(25, 50, 100, 200, 400, 800, 1600)
simrunstotal <- 10000
bootruns <- 500
alphavec <- c(0.01, 0.05, 0.10)
disttable <- data.frame(name = c("normal", "nonnormal"))
disttable$skewness <- list(NULL, c(0.7, 0.9, 1.2, 1.5, 0.8, 1.3))
disttable$kurtosis <- list(NULL, c(3.5, 3.6, 3.7, 3.5, 3.6, 3.7))
######################### Monte Carlo Simulation ###############################
cl <- parallel::makeCluster(nkernel, outfile = "errorcluster.txt")
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/jab49wd/R-projects/R/library"))

simresults <- foreach(jj = 1:nrow(simModels), .packages = c("lavaan", "foreach", "dplyr"), .combine = "rbind", .export = c("jacknife")) %:%
  foreach(n = nobs, .combine = "rbind") %:%
  foreach(sim_runs = 1:simrunstotal, .combine = "rbind") %dopar%
  {
    seed <- round(runif(1, min = 0, max = 100000)*1000, digits = 0)
    correlation <- simModels$correlation[jj]
    temp <- foreach(distn = 1:nrow(disttable), .packages = c("lavaan", "foreach", "dplyr"), .combine = "rbind") %do%
      {
        data <- lavaan::simulateData(model = simModels$model[jj],
                                     sample.nobs = n, # Number of observations.
                                     skewness = disttable$skewness[[distn]],
                                     kurtosis = disttable$kurtosis[[distn]],
                                     seed = seed, # Set random seed.
                                     empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
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
                              , coveragecorr = lowerbound < correlation & correlation <  upperbound
                              , coverageone =  lowerbound < 1 & 1 < upperbound
                              , time = time
                              , seed = seed
            )
            row
          }
        res
      }
    temp
  }
closeAllConnections()
write.csv2(simresults, "HTMTsimresults.csv")

