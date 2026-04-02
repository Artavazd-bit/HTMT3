library(covsim)
library(lavaan)
library(dplyr)
library(foreach)

source("arraysimulation/setup.R")

model <- simModels$model[1]

model_df <- lavaanify(model)


data <- lavaan::simulateData(model = model,
                             empirical = TRUE)

set.seed(seed)

data <- lavaan::simulateData(model = simModels$model[5],
                             sample.nobs = 12,
                             seed = NULL,
                             skewness = NULL, 
                             kurtosis = NULL, 
                             empirical = FALSE,
                             return.type = "data.frame"
)

popcov <- cov(data)
calchtmt(data, 
         model = model_est,
         latent1 = "xi_1",
         latent2 = "xi_2",
         scale = FALSE,
         htmt2 = FALSE)


marginsnorm <- lapply(X = sqrt(diag(popcov)), function(X) list(distr = "norm", sd = X))

marginsxi <- lapply(X = c(4, 5, 6, 4, 5, 6), function(X) list(distr = "chisq", df = X))


data <- vita(marginsnorm, sigma.target = popcov, Nmax = 1000)

vine <- vita(marginsxi, sigma.target = popcov, Nmax = 10000)

round(cov(rvinecopulib::rvine(10^5, data))-popcov, 3)


data <- rvinecopulib::rvine(1000, vine)

colnames(data) <- c("x11", "x12", "x13", "x21", "x22", "x23")
cov(data)

calchtmt(data = data, 
         model = model_est, 
         latent1 = "xi_1",
         latent2 = "xi_2",
         scale = FALSE,
         htmt2 = FALSE)

derivhtmt(data = data, 
          model = model_est,
          latent1 = "xi_1",
          latent2 = "xi_2",
          scale = FALSE,
          htmt2 = FALSE)

calcovcov(data)

deltamethod(data, 
            model = model_est,
            alpha = 0.1,
            latent1 = "xi_1",
            latent2 = "xi_2",
            scale = FALSE,
            htmt2 = FALSE)

wrapper(data = data,
        model = model_est,
        alpha = 0.1,
        latent1 = "xi_1",
        latent2 = "xi_2",
        scale = FALSE,
        htmt2 = FALSE,
        nboot = 500)

n_id <- match(n, unique(conditions$n))



seed
set.seed(seed)
data <- rvinecopulib::rvine(n = 50,
                            vine = vine)
