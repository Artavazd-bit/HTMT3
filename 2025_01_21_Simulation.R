library(lavaan)
library(semTools)
library(cSEM)
library(stringr)
library(doParallel)
library(foreach)
library(dplyr)
library(boot)

source("2024_01_08_functions.R")
source("2024_01_10_gradient_analytically_of_htmt.R")
source("2025_19_02_models_3_4.R")

model_est <- '
              #  latent variables
                xi_1 =~ x11 + x12 + x13 + x14
                xi_2 =~ x21 + x22 + x23 
                
                xi_1 ~~ xi_2
              ' 

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

sim_overview <- foreach(jj = 1:nrow(simModels), .packages = c("lavaan", "semTools", "stringr", "boot"), .combine = "rbind") %:%
                foreach(n = c(50, 100, 200, 500), .combine = "rbind") %:%
                foreach(sim_runs = 1:1000, .combine = "rbind") %dopar%
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
                                                   skewness = NULL, # Numeric vector. The skewness values for the observed variables. Defaults to zero.
                                                   kurtosis = NULL, # Numeric vector. The kurtosis values for the observed variables. Defaults to zero.
                                                   seed = seed, # Set random seed.
                                                   empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                                                   return.type = "data.frame",
                                                   return.fit = FALSE, # If TRUE, return the fitted model that has been used to generate the data as an attribute (called "fit"); this may be useful for inspection.
                                                   debug = FALSE, # If TRUE, debugging information is displayed.
                                                   standardized = FALSE # If TRUE, the residual variances of the observed variables are set in such a way such that the model implied variances are unity. This allows regression coefficients and factor loadings (involving observed variables) to be specified in a standardized metric.
                                                   )
                  #Delta 
                  start_time_delta <- Sys.time()
                  vc_r <- calculate_cov_cov(data = data_cfa)
                  gradient_htmt_1 <- calc_grad_htmt_ana(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = FALSE)
                  Gradient_htmt <- as.matrix(gradient_htmt_1$output$gradient)
                  se_htmt_1 = sqrt(t(Gradient_htmt) %*% vc_r %*% Gradient_htmt / n)
                  z_value_htmt_1 = (gradient_htmt_1$HTMT - 1)/se_htmt_1
                  end_time_delta <- Sys.time()
                  
                  delta_delta <- end_time_delta - start_time_delta
                  alpha_vector <- c(0.9, 0.95, 0.99)
                  
                  CI_HTMT_upper = gradient_htmt_1$HTMT + qnorm(p = alpha_vector, mean = 1, sd = se_htmt_1)
                  names(CI_HTMT_upper) <- c("0.9", "0.95", "0.99")
                  
                  #Bootstrapping
                  set.seed(seed)
                  start_time_boot <- Sys.time()
                  bootstrap <- boot(data_cfa, function(data, indices){calc_htmt(data = data[indices,], model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = TRUE, htmt2 = FALSE)}, R = 500)
                  bootstrap_htmt_1_se <- sd(bootstrap$t, na.rm = TRUE)
                  boot_upper_bound <- quantile(bootstrap$t, na.rm = TRUE, probs = alpha_vector)
                  names(boot_upper_bound) <- c("0.9", "0.95", "0.99")
                  end_time_boot <- Sys.time()
                  
                  delta_boot <- end_time_boot - start_time_boot
                  
                  bootstrap_htmt_1_bias <- bootstrap$t0 - mean(bootstrap$t, na.rm = TRUE)
                  
                  z_value_htmt_1_bootstrap <- (gradient_htmt_1$HTMT - 1) / bootstrap_htmt_1_se
                  
                  save <- data.frame( loading1 = simModels$loading_1[jj]
                                      , loading2 = simModels$loading_2[jj]
                                      , correlation = simModels$correlation[jj]
                                      , n = n
                                      , sim_runs
                                      , htmt_1 = gradient_htmt_1$HTMT
                                      , se_htmt_1_delta = se_htmt_1
                                      , z_value_htmt_1 = z_value_htmt_1
                                      
                                      
                                      , z_test_htmt_09 = z_value_htmt_1 <  qnorm(p = 1 - alpha_vector[1], mean = 0, sd = 1)
                                      , z_test_htmt_095 = z_value_htmt_1 <  qnorm(p = 1 - alpha_vector[2], mean = 0, sd = 1)
                                      , z_test_htmt_099 = z_value_htmt_1 <  qnorm(p = 1 - alpha_vector[3], mean = 0, sd = 1)
                                      
                                      , CI_HTMT_09 = unname(CI_HTMT_upper["0.9"])
                                      , CI_HTMT_095 = unname(CI_HTMT_upper["0.95"])
                                      , CI_HTMT_099 = unname(CI_HTMT_upper["0.99"])
                                      
                                      , CI_HTMT_09_test = unname(CI_HTMT_upper["0.9"]) < 1
                                      , CI_HTMT_095_test = unname(CI_HTMT_upper["0.95"]) < 1
                                      , CI_HTMT_099_test = unname(CI_HTMT_upper["0.99"]) < 1
                      
                                      , se_htmt_1_boot =  bootstrap_htmt_1_se
                                      , boot_htmt_1_bias = bootstrap_htmt_1_bias
                                      , z_value_htmt_1_boot = z_value_htmt_1_bootstrap
                                      
                                      
                                      , z_test_htmt_boot_09 = z_value_htmt_1_bootstrap <  qnorm(p = 1 - alpha_vector[1], mean = 0, sd = 1)
                                      , z_test_htmt_boot_095 = z_value_htmt_1_bootstrap <  qnorm(p = 1 - alpha_vector[2], mean = 0, sd = 1)
                                      , z_test_htmt_boot_099 = z_value_htmt_1_bootstrap <  qnorm(p = 1 - alpha_vector[3], mean = 0, sd = 1)
                                      
                                      , boot_ub_09 = unname(boot_upper_bound["0.9"])
                                      , boot_ub_095 = unname(boot_upper_bound["0.95"])
                                      , boot_ub_099 = unname(boot_upper_bound["0.99"])
                                      
                                      , CI_boot_09 = unname(boot_upper_bound["0.9"]) < 1 
                                      , CI_boot_095 = unname(boot_upper_bound["0.95"]) < 1 
                                      , CI_boot_099 = unname(boot_upper_bound["0.99"]) < 1 
                                      , seed = seed
                                      
                                      , comp_time_delta = delta_delta
                                      , comp_time_boot = delta_boot
                                      )
                  save
                }

closeAllConnections()
sim_overview_without_NA <- sim_overview[complete.cases(sim_overview[,"se_htmt_2"]),]

sim_overview_without_NA <- na.omit(sim_overview)

sim_overview_2 <- sim_overview_without_NA %>% 
  group_by(loading1, loading2, correlation, n) %>%
  summarize(Rejection_rate_htmt_09= mean(z_test_htmt_09), 
            Rejection_rate_htmt_095= mean(z_test_htmt_095),
            Rejection_rate_htmt_099= mean(z_test_htmt_099),
            
            Rejection_CI_09 = mean(CI_HTMT_09_test), 
            Rejection_CI_095 = mean(CI_HTMT_095_test),
            Rejection_CI_099 = mean(CI_HTMT_099_test),
            
            Rejection_boot_09 = mean(CI_boot_09), 
            Rejection_boot_095 = mean(CI_boot_095),
            Rejection_boot_099 = mean(CI_boot_099), 
            
            comp_time_delta = mean(comp_time_delta), 
            comp_time_boot = mean(comp_time_boot)
            )

write.csv2(sim_overview, file= "2025_03_11_sim_overview.csv")
write.csv2(sim_overview_2, file = "2025_03_11_sim_overview2.csv")
