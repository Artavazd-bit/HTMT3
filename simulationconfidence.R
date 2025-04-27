.libPaths(c("/home/jab49wd/R-projects/R/library", .libPaths()))

library(lavaan)
library(doParallel)
library(foreach)
library(dplyr)
library(boot)

source("functionsandmodelsconfidence.R")

##################### Einstellungen ############################################
simModels <- simModels
# nkernel <- 32
nkernel <- 16
# nobs <- c(12, 25, 50, 100, 200, 400, 800)
nobs <- 50
simrunnumer <- 10
bootruns <- 500
######################### Normal data ##########################################
cl <- parallel::makeCluster(nkernel)
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/jab49wd/R-projects/R/library"))

simnormal <- foreach(jj = 1:nrow(simModels), .packages = c("lavaan", "semTools", "stringr", "boot"), .combine = "rbind") %:%
  foreach(n = nobs, .combine = "rbind") %:%
  foreach(sim_runs = 1:simrunnumer, .combine = "rbind") %dopar%
  {
    seed <- round(runif(1, min = 0, max = 100000)*1000, digits = 0)
    data_cfa <- lavaan::simulateData(model = simModels$model[jj],
                                     model.type = "cfa",
                                     meanstructure = FALSE, # means of observed variables enter the model
                                     int.ov.free = FALSE, # if false, intercepts of observed are fixed to zero
                                     int.lv.free = FALSE, # if false, intercepts of latent var fixed to zero
                                     marker.int.zero = FALSE, # only relevant, if the metric of each latent var is set by fixing the first factor loading to unity
                                     conditional.x = FALSE, # If TRUE, we set up the model on the exogenous "x" covariates, the model implied sample statistics only include the non-x variables. If FALSE x are modelled jointly with the other variables and the model implied statistics reflect both sets of variables. 
                                     fixed.x = FALSE, # if TRUE, ex x are considered fixed
                                     orthogonal = FALSE, # if TRUE exogenous latent variables are assumed to be uncorrelated
                                     std.lv = FALSE, # If TRUE, the metric of each latent variable is determined by fixing their variances to 1.0. If FALSE, the metric of each latent variable is determined by fixing the factor loading of the first indicator to 1.0.
                                     auto.fix.first = FALSE, # If TRUE, the factor loading of the first indicator is set to 1.0 for every latent variable
                                     auto.fix.single = FALSE, # If TRUE, the residual variance (if included) of an observed indicator is set to zero if it is the only indicator of a latent variable.
                                     auto.var = TRUE, # If TRUE, the (residual) variances of both observed and latent variables are set free.
                                     auto.cov.lv.x = TRUE, # If TRUE, the covariances of exogenous latent variables are included in the model and set free.
                                     auto.cov.y = FALSE,# If TRUE, the covariances of dependent variables (both observed and latent) are included in the model and set free.
                                     sample.nobs = n, # Number of observations.
                                     ov.var = NULL,# The user-specified variances of the observed variables.
                                     group.label = paste("G", 1:ngroups, sep = ""), # The group labels that should be used if multiple groups are created.
                                     skewness = NULL,
                                     kurtosis = NULL,
                                     seed = seed, # Set random seed.
                                     empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                                     return.type = "data.frame",
                                     return.fit = FALSE, # If TRUE, return the fitted model that has been used to generate the data as an attribute (called "fit"); this may be useful for inspection.
                                     debug = FALSE, # If TRUE, debugging information is displayed.
                                     standardized = FALSE # If TRUE, the residual variances of the observed variables are set in such a way such that the model implied variances are unity. This allows regression coefficients and factor loadings (involving observed variables) to be specified in a standardized metric.
    )
    
    sim <- simFun(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", alpha = c(0.01, 0.05, 0.1), scale = FALSE, htmt2 = FALSE, seed = seed, bootruns = bootruns)
    
    
    delta001 <- data.frame(correlation = simModels$correlation[jj] 
                        , n = n
                        , sim_runs
                        , datatype = "normal"
                        , method = "delta"
                        , alpha = 0.01
                        , lowerbound = sim$delta$lowerbound[1]
                        , upperbound = sim$delta$upperbound[1]
                        , coveragecorr = sim$delta$lowerbound[1] < simModels$correlation[jj] & simModels$correlation[jj] < sim$delta$upperbound[1]
                        , coverageone = sim$delta$lowerbound[1] < 1 & 1 < sim$delta$upperbound[1]
                        , time = sim$delta$time
                        )
    delta005 <- data.frame(correlation = simModels$correlation[jj] 
                       , n = n
                       , sim_runs
                       , datatype = "normal"
                       , method = "delta"
                       , alpha = 0.05
                       , lowerbound = sim$delta$lowerbound[2]
                       , upperbound = sim$delta$upperbound[2]
                       , coveragecorr = sim$delta$lowerbound[2] < simModels$correlation[jj] & simModels$correlation[jj] < sim$delta$upperbound[2]
                       , coverageone = sim$delta$lowerbound[2] < 1 & 1 < sim$delta$upperbound[2]
                       , time = sim$delta$time
                       )
    delta010 <- data.frame(correlation = simModels$correlation[jj] 
                           , n = n
                           , sim_runs
                           , datatype = "normal"
                           , method = "delta"
                           , alpha = 0.10
                           , lowerbound = sim$delta$lowerbound[3]
                           , upperbound = sim$delta$upperbound[3]
                           , coveragecorr = sim$delta$lowerbound[3] < simModels$correlation[jj] & simModels$correlation[jj] < sim$delta$upperbound[3]
                           , coverageone = sim$delta$lowerbound[3] < 1 & 1 < sim$delta$upperbound[3]
                           , time = sim$delta$time
                           )
    
    boot001 <- data.frame(correlation = simModels$correlation[jj] 
                           , n = n
                           , sim_runs
                           , datatype = "normal"
                           , method = "boot"
                           , alpha = 0.01
                           , lowerbound = unname(sim$boot$lowerbound[1])
                           , upperbound = unname(sim$boot$upperbound[1])
                           , coveragecorr = sim$boot$lowerbound[1] < simModels$correlation[jj] & simModels$correlation[jj] < sim$boot$upperbound[1]
                           , coverageone = sim$boot$lowerbound[1] < 1 & 1 < sim$boot$upperbound[1]
                           , time = sim$boot$time
    )
    boot005 <- data.frame(correlation = simModels$correlation[jj] 
                           , n = n
                           , sim_runs
                           , datatype = "normal"
                           , method = "boot"
                           , alpha = 0.05
                           , lowerbound = unname(sim$boot$lowerbound[2])
                           , upperbound = unname(sim$boot$upperbound[2])
                           , coveragecorr = sim$boot$lowerbound[2] < simModels$correlation[jj] & simModels$correlation[jj] < sim$boot$upperbound[2]
                           , coverageone = sim$boot$lowerbound[2] < 1 & 1 < sim$boot$upperbound[2]
                           , time = sim$boot$time
    )
    boot010 <- data.frame(correlation = simModels$correlation[jj] 
                           , n = n
                           , sim_runs
                           , datatype = "normal"
                           , method = "boot"
                           , alpha = 0.10
                           , lowerbound = unname(sim$boot$lowerbound[3])
                           , upperbound = unname(sim$boot$upperbound[3])
                           , coveragecorr = sim$boot$lowerbound[3] < simModels$correlation[jj] & simModels$correlation[jj] < sim$boot$upperbound[3]
                           , coverageone = sim$boot$lowerbound[3] < 1 & 1 < sim$boot$upperbound[3]
                           , time = sim$boot$time
    )
    rbind(delta001
          , delta005
          , delta010
          , boot001
          , boot005
          , boot010
    )
  }

closeAllConnections()


######################### non normal data #####################################
cl <- parallel::makeCluster(nkernel)
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/jab49wd/R-projects/R/library"))

simnonnormal <- foreach(jj = 1:nrow(simModels), .packages = c("lavaan", "semTools", "stringr", "boot"), .combine = "rbind") %:%
  foreach(n = nobs, .combine = "rbind") %:%
  foreach(sim_runs = 1:simrunnumer, .combine = "rbind") %dopar%
  {
    seed <- round(runif(1, min = 0, max = 100000)*1000, digits = 0)
    data_cfa <- lavaan::simulateData(model = simModels$model[jj],
                                     model.type = "cfa",
                                     meanstructure = FALSE, # means of observed variables enter the model
                                     int.ov.free = FALSE, # if false, intercepts of observed are fixed to zero
                                     int.lv.free = FALSE, # if false, intercepts of latent var fixed to zero
                                     marker.int.zero = FALSE, # only relevant, if the metric of each latent var is set by fixing the first factor loading to unity
                                     conditional.x = FALSE, # If TRUE, we set up the model on the exogenous "x" covariates, the model implied sample statistics only include the non-x variables. If FALSE x are modelled jointly with the other variables and the model implied statistics reflect both sets of variables. 
                                     fixed.x = FALSE, # if TRUE, ex x are considered fixed
                                     orthogonal = FALSE, # if TRUE exogenous latent variables are assumed to be uncorrelated
                                     std.lv = FALSE, # If TRUE, the metric of each latent variable is determined by fixing their variances to 1.0. If FALSE, the metric of each latent variable is determined by fixing the factor loading of the first indicator to 1.0.
                                     auto.fix.first = FALSE, # If TRUE, the factor loading of the first indicator is set to 1.0 for every latent variable
                                     auto.fix.single = FALSE, # If TRUE, the residual variance (if included) of an observed indicator is set to zero if it is the only indicator of a latent variable.
                                     auto.var = TRUE, # If TRUE, the (residual) variances of both observed and latent variables are set free.
                                     auto.cov.lv.x = TRUE, # If TRUE, the covariances of exogenous latent variables are included in the model and set free.
                                     auto.cov.y = FALSE,# If TRUE, the covariances of dependent variables (both observed and latent) are included in the model and set free.
                                     sample.nobs = n, # Number of observations.
                                     ov.var = NULL,# The user-specified variances of the observed variables.
                                     group.label = paste("G", 1:ngroups, sep = ""), # The group labels that should be used if multiple groups are created.
                                     skewness = c(0.7, 0.9, 1.2, 1.5, 0.8, 1.3),
                                     kurtosis = c(3.5, 3.6, 3.7, 3.5, 3.6, 3.7),
                                     seed = seed, # Set random seed.
                                     empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                                     return.type = "data.frame",
                                     return.fit = FALSE, # If TRUE, return the fitted model that has been used to generate the data as an attribute (called "fit"); this may be useful for inspection.
                                     debug = FALSE, # If TRUE, debugging information is displayed.
                                     standardized = FALSE # If TRUE, the residual variances of the observed variables are set in such a way such that the model implied variances are unity. This allows regression coefficients and factor loadings (involving observed variables) to be specified in a standardized metric.
    )
    
    sim <- simFun(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", alpha = c(0.01, 0.05, 0.1), scale = FALSE, htmt2 = FALSE, seed = seed, bootruns = bootruns)
    
    
    delta001 <- data.frame(correlation = simModels$correlation[jj] 
                           , n = n
                           , sim_runs
                           , datatype = "nonnormal"
                           , method = "delta"
                           , alpha = 0.01
                           , lowerbound = sim$delta$lowerbound[1]
                           , upperbound = sim$delta$upperbound[1]
                           , coveragecorr = sim$delta$lowerbound[1] < simModels$correlation[jj] & simModels$correlation[jj] < sim$delta$upperbound[1]
                           , coverageone = sim$delta$lowerbound[1] < 1 & 1 < sim$delta$upperbound[1]
                           , time = sim$delta$time
    )
    delta005 <- data.frame(correlation = simModels$correlation[jj] 
                           , n = n
                           , sim_runs
                           , datatype = "nonnormal"
                           , method = "delta"
                           , alpha = 0.05
                           , lowerbound = sim$delta$lowerbound[2]
                           , upperbound = sim$delta$upperbound[2]
                           , coveragecorr = sim$delta$lowerbound[2] < simModels$correlation[jj] & simModels$correlation[jj] < sim$delta$upperbound[2]
                           , coverageone = sim$delta$lowerbound[2] < 1 & 1 < sim$delta$upperbound[2]
                           , time = sim$delta$time
    )
    delta010 <- data.frame(correlation = simModels$correlation[jj] 
                           , n = n
                           , sim_runs
                           , datatype = "nonnormal"
                           , method = "delta"
                           , alpha = 0.10
                           , lowerbound = sim$delta$lowerbound[3]
                           , upperbound = sim$delta$upperbound[3]
                           , coveragecorr = sim$delta$lowerbound[3] < simModels$correlation[jj] & simModels$correlation[jj] < sim$delta$upperbound[3]
                           , coverageone = sim$delta$lowerbound[3] < 1 & 1 < sim$delta$upperbound[3]
                           , time = sim$delta$time
    )
    
    boot001 <- data.frame(correlation = simModels$correlation[jj] 
                          , n = n
                          , sim_runs
                          , datatype = "nonnormal"
                          , method = "boot"
                          , alpha = 0.01
                          , lowerbound = unname(sim$boot$lowerbound[1])
                          , upperbound = unname(sim$boot$upperbound[1])
                          , coveragecorr = sim$boot$lowerbound[1] < simModels$correlation[jj] & simModels$correlation[jj] < sim$boot$upperbound[1]
                          , coverageone = sim$boot$lowerbound[1] < 1 & 1 < sim$boot$upperbound[1]
                          , time = sim$boot$time
    )
    boot005 <- data.frame(correlation = simModels$correlation[jj] 
                          , n = n
                          , sim_runs
                          , datatype = "nonnormal"
                          , method = "boot"
                          , alpha = 0.05
                          , lowerbound = unname(sim$boot$lowerbound[2])
                          , upperbound = unname(sim$boot$upperbound[2])
                          , coveragecorr = sim$boot$lowerbound[2] < simModels$correlation[jj] & simModels$correlation[jj] < sim$boot$upperbound[2]
                          , coverageone = sim$boot$lowerbound[2] < 1 & 1 < sim$boot$upperbound[2]
                          , time = sim$boot$time
    )
    boot010 <- data.frame(correlation = simModels$correlation[jj] 
                          , n = n
                          , sim_runs
                          , datatype = "nonnormal"
                          , method = "boot"
                          , alpha = 0.10
                          , lowerbound = unname(sim$boot$lowerbound[3])
                          , upperbound = unname(sim$boot$upperbound[3])
                          , coveragecorr = sim$boot$lowerbound[3] < simModels$correlation[jj] & simModels$correlation[jj] < sim$boot$upperbound[3]
                          , coverageone = sim$boot$lowerbound[3] < 1 & 1 < sim$boot$upperbound[3]
                          , time = sim$boot$time
    )
    rbind(delta001
          , delta005
          , delta010
          , boot001
          , boot005
          , boot010
    )
  }

closeAllConnections()



simresults <- rbind(simnormal, simnonnormal)

write.csv2(simresults, "simresultsconfidence.csv")






