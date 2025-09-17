library(lavaan)
library(githubinstall)
library(numDeriv)

source("2024_01_08_functions.R")
source("2024_01_10_gradient_analytically_of_htmt.R")
source("2025_19_02_models_3_4.R")


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
                                 empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                                 return.type = "data.frame",
                                 return.fit = FALSE, # If TRUE, return the fitted model that has been used to generate the data as an attribute (called "fit"); this may be useful for inspection.
                                 debug = FALSE, # If TRUE, debugging information is displayed.
                                 standardized = FALSE # If TRUE, the residual variances of the observed variables are set in such a way such that the model implied variances are unity. This allows regression coefficients and factor loadings (involving observed variables) to be specified in a standardized metric.
)



head(data_cfa)

# first i want to see how the data fits the model
fit_cfa <- lavaan::cfa(model  = model_est, 
                       data = data_cfa)

datscale <- scale(data_cfa)

fit_cfa2 <- lavaan::cfa(model  = model_est, 
                       data = datscale)

summary(fit_cfa2)

lavdata <- fit_cfa@Data
lavoptions <- lavInspect(fit_cfa, "options")
lavoptions$gamma.unbiased <- FALSE
lavoptions$gamma.n.minus.one <- FALSE
lavoptions$correlation <- FALSE
lavoptions$conditional.x <- FALSE
lavoptions$fixed.x <- FALSE
lavoptions$meanstructure <- FALSE

test <- lavaan::lav_samplestats_from_data(lavdata = lavdata, lavoptions = lavoptions, NACOV = TRUE)

#ahhh er nimmt die varianz mit, deswegen hier eine 28x28 Matrix.
test@NACOV

cov(data_cfa)*((100-1)/100) - test@cov[[1]]
# ergibt gleich 0

test@NACOV[[1]]

out <- deltamethod(data = data_cfa, model = model_est, alpha = 0.05, latent1 = "xi_1", latent2 = "xi_2", scale = FALSE, htmt2 = FALSE)
vc_r[,1] 

identical(out$omega, test@NACOV[[1]])

matrix_A <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 3, ncol = 4)
matrix_B <- matrix(c(9, 8, 7, 6, 5, 4, 3, 2, 1), nrow = 3, ncol = 3)
crossprod(matrix_A, matrix_B)

# das sind einfach meine Daten
nacov <- function(data){
  Y <- unname(as.matrix(data))
  #zentrieren der Daten
  Yc <- t(t(Y) - colMeans(Y, na.rm = TRUE))
  rows <- ncol(Y)
  z <- sequence(rows)
  idx1 = rep(z[-length(z)], times = rev(tail(z, -1))-1)
  idx2 = unlist(lapply(2:rows, function(x) x:rows), use.names = FALSE)
  # berechnung von Z ??????
  Z <- (Yc[, idx1, drop = FALSE] * Yc[, idx2, drop = FALSE])
  # zentrierung der Daten??????
  Zc <- t(t(Z) - colMeans(Z, na.rm = TRUE))
  #das ist das gleiche wie test@NACOV
  NACOV <- base::crossprod(Zc) / nrow(data)
  NACOV
}
all.equal(nacov(data_cfa), out$omega)

# herausfinden der Indize der unteren Dreiecksmatrix


vc_rtest <- calculate_cov_cov(data = data_cfa)

all.equal(zz, vc_rtest)

########### Test gradient ##################
# funktioniert einwandfrei!1!!
# einziges Problem besteht noch bei Varianz Covarianz matrix der covarianzen. 
htmt1 <- function(x){
  cor_values = wis(cov(data_cfa), model = model_est, latent1 = "xi_1", latent = "xi_2")
  cor_values$val = x  
  calchtmt(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2",scale = FALSE, htmt2 = FALSE)
}


wis <- function(covdata, model, latent1, latent2){
  indicators <- extract_indicators(lv1 = latent1, lv2 = latent2, model_syntax = model)
  
  all_indicators <- unlist(indicators)
  
  cor_subset_data <- covdata
  ind <- which( lower.tri(cor_subset_data,diag=F) , arr.ind = TRUE )
  cor_values <- data.frame( col = dimnames(cor_subset_data)[[2]][ind[,2]] ,
                            row = dimnames(cor_subset_data)[[1]][ind[,1]] ,
                            val = cor_subset_data[ ind ] )
  
  cor_values$type[cor_values$col %in% unlist(indicators[1]) & cor_values$row %in% unlist(indicators[1])] <- "mono1"
  cor_values$type[cor_values$col %in% unlist(indicators[2]) & cor_values$row %in% unlist(indicators[2])] <- "mono2"
  cor_values$type[cor_values$col %in% unlist(indicators[1]) & cor_values$row %in% unlist(indicators[2])] <- "het"
  
  cor_values
}

y <- wis(cov(data_cfa), model = model_est, latent1 = "xi_1", latent = "xi_2")
grad <- numDeriv::grad(func = htmt1, x = y$val) 

calc_grad_htmt_ana(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = FALSE)

# nicht vergessen, semTools berechnet den HTMT auf Cor Basis nicht Cov Basis. selbst wenn man cov() übergibt????
htmt <- semTools::htmt(model = model_est,
                       sample.cov = cov(data_cfa),
                       data = NULL,
                       htmt2 = FALSE,
                       absolute = FALSE
)

calc_htmt(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = FALSE, htmt2 = FALSE)
