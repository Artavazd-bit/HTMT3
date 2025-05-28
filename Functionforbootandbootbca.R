################################################################################
##                                BootBcA                                     
################################################################################
source("functionsandmodelsconfidence.R")
library(dplyr)
library(foreach)
library(doParallel)
################################################################################
## Normal Bootstrap
################################################################################
data <- readRDS("exampledata.rds")
bootruns <- 10

HTMT_boot <- data.frame(matrix(ncol = 2, nrow = 0))
name <- c("HTMT_val", "Boot_seed")
colnames(HTMT_boot) <- name

for(i in 1:bootruns){
  seed <- round(runif(1, min = 0, max = 100000)*1000, digits = 0)
  set.seed(seed)
  datanew <- sample_n(data, nrow(data), replace = TRUE,)
  HTMT_boot[i, "HTMT_val"] <- calchtmt(data = datanew, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = FALSE, htmt2 = FALSE)
  HTMT_boot[i, "Boot_seed"] <- seed
  rm(datanew, seed, i)
}

statistic <- calchtmt(data = data, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = FALSE, htmt2 = FALSE)
HTMT_boot$check <- HTMT_boot$HTMT_val < HTMToriginal
z0 <- qnorm(p = mean(HTMT_boot$check), mean = 0, sd = 1)
################################################################################
## Jackknife
################################################################################
HTMT_jackknife <- data.frame(matrix(ncol = 1, nrow = 0))
colnames(HTMT_jackknife) <- c("HTMT_val")
for(i in 1:nrow(data)){
  datanew <- data[-i,]
  HTMT_jackknife[i, "HTMT_val"] <- calchtmt(data = datanew, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = FALSE, htmt2 = FALSE)
  rm(datanew, i)
}

meanHTMTjacked <- mean(HTMT_jackknife[, "HTMT_val"])
HTMT_jackknife$centered <- HTMT_jackknife$HTMT_val - meanHTMTjacked
HTMT_jackknife$centeredtri <- HTMT_jackknife$centered^3
HTMT_jackknife$centereddos <- HTMT_jackknife$centered^2

accelerator <- sum(HTMT_jackknife$centeredtri) / (6* (sum(HTMT_jackknife$centereddos))^(3/2))

################################################################################
## bca
################################################################################


################################################################################
##              Funktion boot
################################################################################
bootfun <- function(data, nboot, alpha = 0.05, statisticfun, ..., parallel = TRUE, cores = NULL){
  if (parallel) 
    {
    if (is.null(cores)) 
      {
      cores <- parallel::detectCores() - 1
      }
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }
  
  
  if(parallel)
  {
  bootall <- foreach(i = 1:nboot, .combine = "rbind", .packages = c("dplyr", "lavaan", "semTools", "stringr")) %dopar%
    {
      seed <- round(runif(1, min = 0, max = 100000)*1000, digits = 0)
      set.seed(seed)
      datanew <- dplyr::sample_n(data, nrow(data), replace = TRUE)
      tryCatch({
        statistic <- statisticfun(datanew, ...)
        data.frame(
          iteration = i,
          stat = statistic,
          seed = seed
        )
      }, error = function(e){
        warning(paste("Fehler in Iteration", i , ":", e$message))
        data.frame(
            iteration = i,
            stat = NA,
            seed = seed
        )
      })
    }
  }else{
    bootall <- foreach(i = 1:nboot, .combine = "rbind", .packages = c("dplyr", "lavaan", "stringr")) %do%
      {
        seed <- round(runif(1, min = 0, max = 100000)*1000, digits = 0)
        set.seed(seed)
        datanew <- dplyr::sample_n(data, nrow(data), replace = TRUE)
        tryCatch({
          statistic <- statisticfun(datanew, ...)
          data.frame(
            iteration = i,
            stat = statistic,
            seed = seed
          )
        }, error = function(e){
          warning(paste("Fehler in Iteration", i , ":", e$message))
          data.frame(
            iteration = i,
            stat = NA,
            seed = seed
          )
        })
      }
  }
  
  if (parallel) {
    parallel::stopCluster(cl)
    closeAllConnections()
  }
  
  
  
  valid_stats <- bootall$stat[!is.na(bootall$stat)]
  
  lowerboundb <- unname(quantile(valid_stats, probs = alpha/2))
  upperboundb <- unname(quantile(valid_stats, probs = 1 - (alpha/2)))
  
  bootmean <- mean(valid_stats)
  bootsd <- sd(valid_stats)
  
  return(list(
    table = bootall,
    lowerbound = lowerboundb, 
    upperbound = upperboundb,
    mean = bootmean, 
    sd = bootsd,
    n_valid = length(valid_stats), 
    n_total = nboot,
    alpha = alpha
  ))
}

boot <- bootfun(data = data, nboot = 10, alpha = 0.05, statisticfun = calchtmt, model = model_est, latent1 = "xi_1", 
        latent2 <- "xi_2", scale = FALSE, htmt2 = FALSE, parallel = TRUE)
################################################################################
## jacknife
################################################################################
jacknifefun <- function(data, statisticfun, ..., parallel = TRUE, cores = NULL, alpha = 0.05)
{
  if (parallel) 
  {
    if (is.null(cores)) 
    {
      cores <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }
  
  
  if(parallel)
  {
    jacknifetab <- foreach(i = 1:nrow(data), .combine = "rbind", .packages = c("dplyr", "lavaan", "stringr")) %dopar%
    {
      datanew <- data[-i,]
      tryCatch({
        statistic <- statisticfun(datanew, ...)
        data.frame(
          iteration = i,
          stat = statistic
        )
      }, error = function(e){
        warning(paste("Fehler in Iteration", i , ":", e$message))
        data.frame(
          iteration = i,
          stat = NA
        )
      })
    }
  }else{
    jacknifetab <- foreach(i = 1:nrow(data), .combine = "rbind", .packages = c("dplyr", "lavaan", "stringr")) %do%
      {
        datanewtab <- data[-i,]
        tryCatch({
          statistic <- statisticfun(datanew, ...)
          data.frame(
            iteration = i,
            stat = statistic
          )
        }, error = function(e){
          warning(paste("Fehler in Iteration", i , ":", e$message))
          data.frame(
            iteration = i,
            stat = NA
          )
        })
      }
  }
  
  if (parallel) {
    parallel::stopCluster(cl)
    closeAllConnections()
  }
  
  valid_stats <- jacknifetab$stat[!is.na(jacknifetab$stat)]
  
  lowerboundb <- unname(quantile(valid_stats, probs = alpha/2))
  upperboundb <- unname(quantile(valid_stats, probs = 1 - (alpha/2)))
  
  jackmean <- mean(valid_stats)
  jacksd <- sd(valid_stats)
  
  statcentered <- jacknifetab$stat - jackmean
  statcenteredsq <- (jacknifetab$stat - jackmean)^2
  statcenteredthree <- (jacknifetab$stat - jackmean)^3
  
  sumthree <- sum(statcenteredthree)
  sumsq <- sum(statcenteredsq)
  
  accelerator <- sumthree / (6*sumsq^(3/2))
  
  return(list(
    table = jacknifetab,
    lowerbound = lowerboundb, 
    upperbound = upperboundb,
    mean = jackmean, 
    sd = jacksd,
    n_valid = length(valid_stats), 
    alpha = alpha,
    accelerator = accelerator
  ))
}
jacknife <- jacknifefun(data = data, statisticfun = calchtmt, model = model_est, latent1 = "xi_1", 
            latent2 <- "xi_2", scale = FALSE, htmt2 = FALSE, alpha = 0.05, parallel = TRUE)
################################################################################
bootbcafun <- function(data, nboot, alpha = 0.05, statisticfun, ..., parallel = TRUE, cores = NULL){
  statistic <- statisticfun(data, ...)
  boot <- bootfun(data, nboot, alpha = alpha, statisticfun, ..., parallel = parallel, cores = cores)
  #print(boot)
  z0 <- qnorm(p = mean(boot$table$stat < statistic), mean = 0, sd = 1)
  
  jacknife <- jacknifefun(data, statisticfun, ..., parallel = parallel, cores = cores, alpha = alpha)
  #print(jacknife)
  acc <- jacknife$accelerator
  
  zalpha <- qnorm(p= alpha)
  zalpha2 <- qnorm(p= 1-alpha)
  
  a1 <- pnorm(q = (z0 + (z0 + zalpha)/(1-acc*(z0+zalpha))), sd = 1, mean = 0)
  a2 <- pnorm(q = (z0 + (z0 + zalpha2)/(1-acc*(z0+zalpha2))), sd = 1, mean = 0)
  
  lowerbound <- unname(quantile(boot$table$stat, probs = a1))
  upperbound <- unname(quantile(boot$table$stat, probs = a2))
  
  return(list(lowerbound = lowerbound, upperbound = upperbound))
} 
bootbcafun(data = data, nboot = 10, alpha = 0.05, statisticfun = calchtmt, model = model_est, latent1 = "xi_1", 
latent2 <- "xi_2", scale = FALSE, htmt2 = FALSE, parallel = TRUE)

