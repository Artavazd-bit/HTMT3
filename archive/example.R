
source("2024_01_08_functions.R")
source("2024_01_10_gradient_analytically_of_htmt.R")
model <- '
                xi_1 =~ 0.7*x11 + 0.7*x12 + 0.7*x13
                xi_2 =~ 0.8*x21 + 0.8*x22 + 0.8*x23 + 0.8*x24
                x11 ~~ 2*x11 + 0*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23 + 0*x24
                x12 ~~ 3*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23 + 0*x24
                x13 ~~ 4*x13 + 0*x21 + 0*x22 + 0*x23 + 0*x24
                x21 ~~ 5*x21 + 0*x22 + 0*x23 + 0*x24
                x22 ~~ 6*x22 + 0*x23 + 0*x24
                x23 ~~ 7*x23 + 0*x24
                x24 ~~ 8*x24
              
                xi_1 ~~ 1*xi_1 + 0.9*xi_2
                xi_2 ~~ 1*xi_2
              ' 

model_est <- '
              #  latent variables
                xi_1 =~ x11 + x12 + x13 
                xi_2 =~ x21 + x22 + x23 + x24
                
                xi_1 ~~ xi_2
              ' 

data_cfa <- lavaan::simulateData(model = simModels$model[40], 
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
                                 sample.nobs = 10000, # Number of observations.
                                 ov.var = NULL,# The user-specified variances of the observed variables.
                                 group.label = paste("G", 1:ngroups, sep = ""), # The group labels that should be used if multiple groups are created.
                                 skewness = NULL, # Numeric vector. The skewness values for the observed variables. Defaults to zero.
                                 kurtosis = NULL, # Numeric vector. The kurtosis values for the observed variables. Defaults to zero.
                                 
                                 seed = NULL, # Set random seed.
                                 
                                 empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                                 
                                 return.type = "data.frame",
                                 
                                 return.fit = FALSE, # If TRUE, return the fitted model that has been used to generate the data as an attribute (called "fit"); this may be useful for inspection.
                                 debug = FALSE, # If TRUE, debugging information is displayed.
                                 standardized = FALSE # If TRUE, the residual variances of the observed variables are set in such a way such that the model implied variances are unity. This allows regression coefficients and factor loadings (involving observed variables) to be specified in a standardized metric.
)

n <- nrow(data_cfa)
calc_htmt(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = TRUE, htmt2 = FALSE)
calc_htmt(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = FALSE, htmt2 = FALSE)

cov_cov <- calculate_cov_cov(data = data_cfa)
htmt_ana_cov <- calc_grad_htmt_ana(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = FALSE)
se_htmt_1 = sqrt(t(htmt_ana_cov$output$gradient) %*% cov_cov %*% htmt_ana_cov$output$gradient / n)
se_htmt_1 


cor_cov <- calculate_corr_cov_fast(data = data_cfa)
htmt_ana_cor <- calc_grad_htmt_ana(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = TRUE)
se_htmt_2 <-  sqrt(t(htmt_ana_cor$output$gradient) %*% cor_cov %*% htmt_ana_cor$output$gradient / n)
se_htmt_2

se_htmt_1 - se_htmt_2
