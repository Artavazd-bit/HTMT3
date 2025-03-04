# first test with a reflective two factor model
# one latent variable will have 3 items the other 4 for the beginning 
# later i want to test different  combinations
# i use lavaan::simulateData to simulate the data
# 
library(lavaan)
library(semTools)
library(cSEM)
library(stringr)

source("2024_01_08_functions.R")
source("2024_01_10_gradient_analytically_of_htmt.R")

model_dgp <- '
              #  latent variables
                xi_1 =~ 0.8*x11 + 0.8*x12 + 0.8*x13
                xi_2 =~ 0.7*x21 + 0.7*x22 + 0.7*x23 
                x11 ~~ 1*x11 + 0*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23 
                x12 ~~ 1*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23
                x13 ~~ 1*x13 + 0*x21 + 0*x22 + 0*x23 
                x21 ~~ 1*x21 + 0*x22 + 0*x23
                x22 ~~ 1*x22 + 0*x23 
                x23 ~~ 1*x23
                
              #  fix covariances between xi_1 and xi_2 und setze die Varianz auf 1
                xi_1 ~~ 1*xi_1 + 1*xi_2
                xi_2 ~~ 1*xi_2
              ' 

model_dgp2 <- '
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



model_est <- '
              #  latent variables
                xi_1 =~ x11 + x12 + x13 + x14
                xi_2 =~ x21 + x22 + x23
                
                xi_1 ~~ xi_2
              ' 




data_cfa <- lavaan::simulateData(model = simModels$model[1], 
                            model.type = "cfa",
                            meanstructure = FALSE, # means of observed variables enter the model
                            int.ov.free = FALSE, # if false, intercepts of observed are fixed to zero
                            int.lv.free = FALSE, # if false, intercepts of latent var fixed to zero
                            marker.int.zero = FALSE, # only relevant, if the metric of each latent var is set by fixing the first factor loading to unity
                            conditional.x = FALSE, # If TRUE, we set up the model on the exogenous "x" covariates, the model implied sample statistics only include the non-x variables. If FALSE x are modelled jointly with the other variables and the model implied statistics reflect both sets of variables. 
                            fixed.x = FALSE, # if TRUE, ex x are considered fixed
                            orthogonal = FALSE, # if TRUE exogenous latent variables are assumed to be uncorrelated
                            std.lv = TRUE, # If TRUE, the metric of each latent variable is determined by fixing their variances to 1.0. If FALSE, the metric of each latent variable is determined by fixing the factor loading of the first indicator to 1.0.
                            auto.fix.first = FALSE, # If TRUE, the factor loading of the first indicator is set to 1.0 for every latent variable
                            auto.fix.single = FALSE, # If TRUE, the residual variance (if included) of an observed indicator is set to zero if it is the only indicator of a latent variable.
                            auto.var = TRUE, # If TRUE, the (residual) variances of both observed and latent variables are set free.
                            auto.cov.lv.x = FALSE, # If TRUE, the covariances of exogenous latent variables are included in the model and set free.
                            auto.cov.y = FALSE,# If TRUE, the covariances of dependent variables (both observed and latent) are included in the model and set free.
                            sample.nobs = 100, # Number of observations.
                            ov.var = NULL,# The user-specified variances of the observed variables.
                            group.label = paste("G", 1:ngroups, sep = ""), # The group labels that should be used if multiple groups are created.
                            skewness = NULL, # Numeric vector. The skewness values for the observed variables. Defaults to zero.
                            kurtosis = NULL, # Numeric vector. The kurtosis values for the observed variables. Defaults to zero.
                            
                            seed = NULL, # Set random seed.
                            
                            empirical = TRUE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                            
                            return.type = "data.frame",
                            
                            return.fit = FALSE, # If TRUE, return the fitted model that has been used to generate the data as an attribute (called "fit"); this may be useful for inspection.
                            debug = FALSE, # If TRUE, debugging information is displayed.
                            standardized = FALSE # If TRUE, the residual variances of the observed variables are set in such a way such that the model implied variances are unity. This allows regression coefficients and factor loadings (involving observed variables) to be specified in a standardized metric.
                            )

# first i want to see how the data fits the model
fit_cfa <- lavaan::cfa(model  = model_est, 
                       data = data_cfa
                       )

# i found this function, still need to figure it out
discriminantValidity(object = fit_cfa,
                     cutoff = 0.9,
                     merge = FALSE,
                     level = 0.95
                     )

# data to covariance
cov_data_cfa <- var(data_cfa)

# calculate the htmt of the data from above
htmt <- semTools::htmt(model = model_est,
               data =  data_cfa, 
               sample.cov = NULL,
               htmt2 = FALSE,
              absolute = FALSE
               )
# calculate the covariance-variance matrix of the correlationmatrix of the data from above
# functions are defined in 2024_08_01_functions.R
vc_r <- calculate_corr_cov_fast(data = data_cfa)

# calculate the gradient of the htmt - Function
# it has two outputs in a list, a gradient for htmt and htmt2 
gradient <- calc_gradient(data = data_cfa, # data
                          sim_runs = 1, # param for simulation
                          n = length(data_cfa), 
                          jj = 1, # param for simulation
                          model_est, # model
                          latent1_index = 1, # when i want to include more latent variables i can specify here which relationship i want to test 
                          latent2_index = 2) # when i want to include more latent variables i can specify here which relationship i want to test 

gradient$htmt2$`1e-05`

gr <- t(gradient$htmt$`1e-05`)

htmt_ana_cor <- calc_grad_htmt_ana(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = TRUE)
htmt2_ana_cor <- calc_grad_htmt2_ana(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = TRUE)

htmt_ana_cov <- calc_grad_htmt_ana(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = FALSE)
htmt2_ana_cov <- calc_grad_htmt2_ana(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = FALSE)

sqrt(diag( gr %*% vc_r %*% t(gr)))

# i want to test if i'm siginificantly from 1 different, no ? 
(htmt[1,2] - 1)/sqrt(gr %*% vc_r %*% t(gr))


cSEM::csem(.data = data_cfa, .model = model_est, )


cSEM::calculateHTMT()
calc_htmt(data = data_cfa, model = model_est, latent1 = "xi_1", "xi_2", scale = TRUE, htmt2 = FALSE)

calc_htmt(data = data_cfa, model = model_est, latent1 = "xi_1", "xi_2", scale = TRUE, htmt2 = TRUE)
