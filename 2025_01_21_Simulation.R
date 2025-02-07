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


model_dgp_1 <- '
              #  latent variables
                xi_1 =~ 0.8*x11 + 0.8*x12 + 0.8*x13
                xi_2 =~ 0.7*x21 + 0.7*x22 + 0.7*x23 + 0.7*x24
                x11 ~~ 1*x11 + 0*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23 + 0*x24
                x12 ~~ 1*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23 + 0*x24
                x13 ~~ 1*x13 + 0*x21 + 0*x22 + 0*x23 + 0*x24
                x21 ~~ 1*x21 + 0*x22 + 0*x23 + 0*x24
                x22 ~~ 1*x22 + 0*x23 + 0*x24
                x23 ~~ 1*x23 + 0*x24
                x24 ~~ 1*x24
              #  fix covariances between xi_1 and xi_2 und setze die Varianz auf 1
                xi_1 ~~ 1*xi_1 + 1*xi_2
                xi_2 ~~ 1*xi_2
              ' 

model_dgp_2 <- '
              #  latent variables
                xi_1 =~ 0.8*x11 + 0.8*x12 + 0.8*x13
                xi_2 =~ 0.7*x21 + 0.7*x22 + 0.7*x23 + 0.7*x24
                x11 ~~ 1*x11 + 0*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23 + 0*x24
                x12 ~~ 1*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23 + 0*x24
                x13 ~~ 1*x13 + 0*x21 + 0*x22 + 0*x23 + 0*x24
                x21 ~~ 1*x21 + 0*x22 + 0*x23 + 0*x24
                x22 ~~ 1*x22 + 0*x23 + 0*x24
                x23 ~~ 1*x23 + 0*x24
                x24 ~~ 1*x24
              #  fix covariances between xi_1 and xi_2
                xi_1 ~~ 1*xi_1 + 0.5*xi_2
                xi_2 ~~ 1*xi_2
              ' 
corr_vector <- c(1, 0.5)
model_list <- list(model_dgp_1, model_dgp_2)

model_est <- '
              #  latent variables
                xi_1 =~ x11 + x12 + x13
                xi_2 =~ x21 + x22 + x23 + x24
                
                xi_1 ~~ xi_2
              ' 

HTMT_function <- function(data, indices){
  d <- data[indices,]
  output <- calc_htmt(data = d, model = model_est, latent1 = "xi_1", "xi_2", scale = TRUE, htmt2 = FALSE)
  return(output)
}

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

sim_overview <- foreach(jj = 1:length(model_list), .packages = c("lavaan", "semTools", "stringr", "boot"), .combine = "rbind") %:%
                foreach(n = c(12, 25, 50, 100, 200,  500, 1000), .combine = "rbind") %:%
                foreach(sim_runs = 1:1000, .combine = "rbind") %dopar%
                {
                  seed <- round(runif(1, min = 0, max = n) * 1000, digits = 0)
                  data_cfa <- lavaan::simulateData(model = model_list[[jj]],
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
                  
                  vc_r <- calculate_corr_cov_fast(data = data_cfa)
                  gradient_htmt_1 <- calc_grad_htmt_ana(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = TRUE)
                  #gradient_htmt_2 <- calc_grad_htmt2_ana(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = TRUE) 
                  
                  Gradient_htmt <- as.matrix(gradient_htmt_1$output$gradient)
                  #Gradient_htmt_2 <- as.matrix(gradient_htmt_2$output$gradient)
                  
                  se_htmt_1 = sqrt(t(Gradient_htmt) %*% vc_r %*% Gradient_htmt / n)
                  #se_htmt_2 = sqrt(t(Gradient_htmt_2) %*% vc_r %*% Gradient_htmt_2 / 100)
                  
                  t_value_htmt_1 = (gradient_htmt_1$HTMT - 1)/se_htmt_1
                  #t_value_htmt_2 = (gradient_htmt_2$HTMT2 - 1)/se_htmt_2
                  
                  #Bootstrapping
                  #bootstrap <- boot(data_cfa, HTMT_function, R = 100, seed = seed)
                  set.seed(6064)
                  bootstrap <- boot(data_cfa, function(data, indices){calc_htmt(data = data[indices,], model = model_est, latent1 = "xi_1", "xi_2", scale = TRUE, htmt2 = FALSE)}, R = 500)
                  bootstrap_htmt_1_se <- sd(na.omit(bootstrap$t))
                  bootstrap_htmt_1_bias <- bootstrap$t0 - mean(na.omit(bootstrap$t))
                  t_value_htmt_1_bootstrap <- (gradient_htmt_1$HTMT - 1) / bootstrap_htmt_1_se
                  
                  save <- data.frame( true_corr = corr_vector[jj],
                                      n = n,
                                      sim_runs,
                                      htmt_1 = gradient_htmt_1$HTMT,
                                      se_htmt_1_delta = se_htmt_1,
                                      t_value_htmt_1 = t_value_htmt_1,
                                      t_test_htmt_1 = t_value_htmt_1 < qnorm(0.05),
                                      
                                      se_htmt_1_boot =  bootstrap_htmt_1_se,
                                      boot_htmt_1_bias = bootstrap_htmt_1_bias, 
                                      t_value_htmt_1_boot = t_value_htmt_1_bootstrap, 
                                      t_test_htmt_1_boot = t_value_htmt_1_bootstrap < qnorm(0.05),
                                      seed = seed
                                      
                                      #htmt_2 = gradient_htmt_2$HTMT2,
                                      #se_htmt_2 = se_htmt_2,
                                      #t_value_htmt_2 = t_value_htmt_2,
                                      #t_test_htmt_2 = t_value_htmt_2 < qnorm(0.05),
                                      )
                  #save$seed <- list(seed)
                  save
                }

closeAllConnections()
sim_overview_without_NA <- sim_overview[complete.cases(sim_overview[,"se_htmt_2"]),]

sim_overview_without_NA <- na.omit(sim_overview)

sim_overview_2 <- sim_overview %>% 
  group_by(true_corr, n) %>%
  summarize(Rejection_rate_htmt_1= mean(t_test_htmt_1), 
            #Rejection_rate_htmt_2 = mean(t_test_htmt_2), 
            Rejection_rate_htmt_1_boot = mean(t_test_htmt_1_boot))
