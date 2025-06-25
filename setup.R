########################## Script for functions ################################
## Date: 25.06.2025
## Info: 
# This script includes all functions needed for the simulations study to assess 
# the statistical properties of the HTMT introduced by Henseler et al.(2015)
###################### Gradient Calculation HTMT and HTMT2 #####################
#data <- readRDS("exampledata.rds")

################################################################################
## derivhtmt: calculates the gradient of the HTMT and HTMT2 (htmt2) and if 
#             covariances or correlations should be used (scale)
################################################################################
derivhtmt <- function(data, model, latent1, latent2, scale, htmt2){
  
  # select only the indicators of the two selected latent variables
  model_df <- lavaanify(model)
  listind1 <- list(model_df$rhs[model_df$lhs == latent1 & model_df$op == "=~"])
  listind2 <- list(model_df$rhs[model_df$lhs == latent2 & model_df$op == "=~"])
  all_indicators <- unlist(list(listind1, listind2)) 
  
  subset_data <- data[, all_indicators]
  
  # for tau-equivalent assumption use scale = FALSE, for parallel assumption you 
  # can use scale = TRUE (original version of HTMT)
  if(scale == FALSE){
    cor_subset_data <- cov(subset_data)
  } else { 
    cor_subset_data <- cor(subset_data)
  }  
  
  # only the lower triangular matrix of the cor/cov matrix is relevant for htmt 
  # and htmt2
  ind <- which( lower.tri(cor_subset_data,diag=F) , arr.ind = TRUE )
  cor_values <- data.frame( col = dimnames(cor_subset_data)[[2]][ind[,2]] ,
                            row = dimnames(cor_subset_data)[[1]][ind[,1]] ,
                            val = cor_subset_data[ ind ] )
  
  # i need to distinguish between the type of corr/cov values
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row 
                  %in% unlist(listind1)] <- "mono1"
  cor_values$type[cor_values$col %in% unlist(listind2) & cor_values$row 
                  %in% unlist(listind2)] <- "mono2"
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row 
                  %in% unlist(listind2)] <- "het"
  
  K_i <- length(unlist(listind1))
  K_j <- length(unlist(listind2))
  
  tryCatch({
  # Caluclation of HTMT and HTMT2 and (analytical) gradient 
  if (htmt2 == FALSE){
    A = 1/(K_i*K_j) * sum(cor_values$val[cor_values$type == "het"])
    B = 2/(K_i*(K_i-1)) *  sum(cor_values$val[cor_values$type == "mono1"]) 
    C = 2/(K_j*(K_j-1)) *  sum(cor_values$val[cor_values$type == "mono2"]) 
    HTMT <- A / ((B*C)^(1/2))
    
    cor_values$gradient[cor_values$type == "het"] <- (1/(K_i*K_j) )/((B*C)^(1/2))
    cor_values$gradient[cor_values$type == "mono1"] <- -HTMT * 1/(K_i*(K_i-1)) * B^-1
    cor_values$gradient[cor_values$type == "mono2"] <- -HTMT * 1/(K_j*(K_j-1)) * C^-1
  }
  else if(htmt2 == TRUE){
    A =  prod(cor_values$val[cor_values$type == "het"])^(1/(K_i*K_j))
    B =  prod(cor_values$val[cor_values$type == "mono1"])^(2/(K_i*(K_i-1))) 
    C =  prod(cor_values$val[cor_values$type == "mono2"])^(2/(K_j*(K_j-1))) 
    HTMT <- A / ((B*C)^(1/2))
    
    cor_values$gradient[cor_values$type == "het"] <- (1/(K_i*K_j)) * 
      prod(cor_values$val[cor_values$type == "het"])^((1/(K_i*K_j))-1) * 
      prod(cor_values$val[cor_values$type == "het"])/cor_values$val[cor_values$type == "het"] * 
      1/(sqrt((B*C)))
    cor_values$gradient[cor_values$type == "mono1"] <- A * 1/2 * (2/(K_i*(K_i-1))) * 
      prod(cor_values$val[cor_values$type == "mono1"])^((2/(K_i*(K_i-1)))-1) * 
      prod(cor_values$val[cor_values$type == "mono1"])/cor_values$val[cor_values$type == "mono1"] * 
      C * (B*C)^(-3/2) * -1
    cor_values$gradient[cor_values$type == "mono2"] <- A * 1/2 * (2/(K_j*(K_j-1))) * 
      prod(cor_values$val[cor_values$type == "mono2"])^((2/(K_j*(K_j-1)))-1) * 
      prod(cor_values$val[cor_values$type == "mono2"])/cor_values$val[cor_values$type == "mono2"] * 
      B * (B*C)^(-3/2) * -1
  }else{
    print("ERROR")
  }
  } , error = function(e){
    cat("An error occurred:", e$message, "\n")
    return(NA)
  }, warning = function(w){
    cat("A warning occurred:", w$message, "\n")
    return(NA)
  })
  list(output = cor_values, HTMT = HTMT)
} 
################################################################################
## calchtmt: calculates the htmt or htmt2 (htmt2) and if correlations or 
#           covariances should be used.
# i use this function for the bootstrap and bootbca methods. 
# in the simulation derivhtmt also gives me the HTMT and HTMT2 estimates. 
################################################################################
calchtmt <- function(data, model, latent1, latent2, scale, htmt2){
  
  model_df <- lavaanify(model)
  
  listind1 <- list(model_df$rhs[model_df$lhs == latent1 & model_df$op == "=~"])
  listind2 <- list(model_df$rhs[model_df$lhs == latent2 & model_df$op == "=~"])
  
  all_indicators <- unlist(list(listind1, listind2)) 
  
  subset_data <- data[, all_indicators]
  if(scale == FALSE){
    cor_subset_data <- cov(subset_data)
  } else { 
    cor_subset_data <- cor(subset_data)
  }  
  ind <- which( lower.tri(cor_subset_data,diag=F) , arr.ind = TRUE )
  cor_values <- data.frame( col = dimnames(cor_subset_data)[[2]][ind[,2]] ,
                            row = dimnames(cor_subset_data)[[1]][ind[,1]] ,
                            val = cor_subset_data[ ind ] )
  
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row 
                  %in% unlist(listind1)] <- "mono1"
  cor_values$type[cor_values$col %in% unlist(listind2) & cor_values$row 
                  %in% unlist(listind2)] <- "mono2"
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row 
                  %in% unlist(listind2)] <- "het"
  
  K_i <- length(unlist(listind1))
  K_j <- length(unlist(listind2))
  
  tryCatch({
  if(htmt2 == FALSE){
    A = 1/(K_i*K_j) * sum(cor_values$val[cor_values$type == "het"])
    B = 2/(K_i*(K_i-1)) *  sum(cor_values$val[cor_values$type == "mono1"]) 
    C = 2/(K_j*(K_j-1)) *  sum(cor_values$val[cor_values$type == "mono2"]) 
    HTMT <- A / ((B*C)^(1/2))
  }
  else if(htmt2 == TRUE){
    A =  prod(cor_values$val[cor_values$type == "het"])^(1/(K_i*K_j))
    B =  prod(cor_values$val[cor_values$type == "mono1"])^(2/(K_i*(K_i-1))) 
    C =  prod(cor_values$val[cor_values$type == "mono2"])^(2/(K_j*(K_j-1))) 
    HTMT <- A / ((B*C)^(1/2))
  }
  else{
    print("ERROR")
  }
    return(HTMT)
  } , error = function(e){
    cat("An error occurred:", e$message, "\n")
    return(NA)
  }, warning = function(w){
    cat("A warning occurred:", w$message, "\n")
    return(NA)
  })
}
################################################################################
## calcovcov: Caluclates the variance covariance matrix of the off diagonal 
#              non-redundant covariances. 
# i use this function within the delta method function to derive omega for cov-
# matrices
################################################################################
calcovcov <- function(data) {
  n <- nrow(data)
  p <- ncol(data)
  size_n <- (p * p - p) / 2
  
  data_centered <- scale(data, center = TRUE, scale = FALSE)
  
  # Create indices for upper triangle
  indices <- which(lower.tri(matrix(0, p, p)), arr.ind = TRUE)
  
  # Initialize result matrix
  vc_r <- matrix(0, nrow = size_n, ncol = size_n)
  
  # Pre-calculate all possible products of centered variables
  products <- array(0, dim = c(n, p, p))
  for(i in 1:p) {
    for(j in i:p) {
      products[, i, j] <- data_centered[, i] * data_centered[, j] 
      products[, j, i] <- products[, i, j]
    }
  }
  
  # Calculate covariances using vectorized operations
  for(idx1 in 1:nrow(indices)) {
    x <- indices[idx1, "row"] 
    y <- indices[idx1, "col"]
    for(idx2 in idx1:nrow(indices)) {
      z <- indices[idx2, "row"] 
      t <- indices[idx2, "col"]
      # Calculate fourth-order moments using pre-computed products
      omega_xyzt <- mean(products[, x, y] * products[, z, t]) - 
        ( mean(products[, x, y]) * mean(products[, z, t]) )
      vc_r[idx1, idx2] <- omega_xyzt
      vc_r[idx2, idx1] <- omega_xyzt  
    }
  }
  return(vc_r)
}
################################################################################
## calcovcor: Caluclates the variance covariance matrix of the off diagonal 
#              non-redundant correlations.
# omega for correlation matrices
################################################################################
calcovcor <- function(data) {
  n <- nrow(data)
  p <- ncol(data)
  size_n <- (p * p - p) / 2
  
  # Pre-calculate means and centered data
  data_means <- colMeans(data)
  data_centered <- scale(data, center = TRUE, scale = FALSE)
  data_sd <- apply(data, 2, sd)
  
  # Create indices for upper triangle
  indices <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  
  # Initialize result matrix
  vc_r <- matrix(0, nrow = size_n, ncol = size_n)
  
  # Pre-calculate all possible products of centered variables
  products <- array(0, dim = c(n, p, p))
  for(i in 1:p) {
    for(j in i:p) {
      products[, i, j] <- data_centered[, i] * data_centered[, j]
      products[, j, i] <- products[, i, j]
    }
  }
  
  # Function to get position in upper triangular matrix
  get_pos <- function(i, j) {
    if(i > j) {
      temp <- i
      i <- j
      j <- temp
    }
    return(p*(i-1) - i*(i-1)/2 + j - i)
  }
  
  # Calculate covariances using vectorized operations
  for(idx1 in 1:nrow(indices)) {
    x <- indices[idx1, 1]
    y <- indices[idx1, 2]
    
    for(idx2 in idx1:nrow(indices)) {
      z <- indices[idx2, 1]
      t <- indices[idx2, 2]
      
      sd <- sd(data[,x]) * sd(data[,y]) * sd(data[,z]) * sd(data[,t])
      
      
      # Calculate fourth-order moments using pre-computed products
      mu_xyzt <- mean(products[, x, y] * products[, z, t])  / sd
      
      mu_xxzt <- mean(products[, x, x] * products[, z, t]) / sd
      mu_yyzt <- mean(products[, y, y] * products[, z, t]) / sd
      mu_xyzz <- mean(products[, x, y] * products[, z, z]) / sd
      mu_xytt <- mean(products[, x, y] * products[, t, t]) / sd
      
      mu_xxtt <- mean(products[, x, x] * products[, t, t]) / sd
      mu_xxzz <- mean(products[, x, x] * products[, z, z]) / sd
      mu_yytt <- mean(products[, y, y] * products[, t, t]) / sd
      mu_yyzz <- mean(products[, y, y] * products[, z, z]) / sd
      
      # Get correlation coefficients
      rxy <- cor(data[, x], data[, y])
      rzt <- cor(data[, z], data[, t])
      
      # Calculate covariance
      cov_val <- (mu_xyzt - 0.5 * rxy * (mu_xxzt + mu_yyzt) - 
                    0.5 * rzt * (mu_xyzz + mu_xytt) + 
                    0.25 * rxy * rzt * (mu_xxzz + mu_xxtt + mu_yyzz + mu_yytt))
      
      pos1 <- get_pos(x, y)
      pos2 <- get_pos(z, t)
      
      vc_r[pos1, pos2] <- cov_val
      vc_r[pos2, pos1] <- cov_val  # Matrix is symmetric
    }
  }
  
  return(vc_r)
}
################################################################################
## deltamethod: calculates the confidence intervals and standard errors with the
#               help of the delta method. The output is the HTMT, the se and the 
#                lowerbound and upperbound.
################################################################################
deltamethod <- function(data, model, alpha, latent1, latent2, scale, htmt2)
{
  starttime <- Sys.time()
  gdf <- derivhtmt(data = data, model = model, latent1 = latent1, latent2 = latent2, scale = scale, htmt2 = htmt2)
  gradient <- as.matrix(gdf$output$gradient)
  if(scale == FALSE){
    omega <- calcovcov(data = data)
  } else if(scale == TRUE){
    omega <- calcovcor(data = data)
  } else print("ERROR")
  se <- sqrt(t(gradient) %*% omega %*% gradient / nrow(data))[1]
  # zvalue <- (gdf$HTMT - test)/se
  # here i want to test whether im in the lowest alpha percent cases. 
  # ztest <- zvalue <  qnorm(p = alpha, mean = 0, sd = 1)
  upperbound <- gdf$HTMT + qnorm(p = 1 - (alpha/2), mean = 0, sd = 1)*se
  lowerbound <- gdf$HTMT - qnorm(p = 1 - (alpha/2), mean = 0, sd = 1)*se
  endtime <- Sys.time()
  tdelta <- endtime - starttime
  
  if (htmt2 == FALSE){
    list(HTMT = gdf$HTMT, se = se, lowerbound = lowerbound, upperbound = upperbound, time = tdelta, missing = NA)
  }else if (htmt2 == TRUE){
    list(HTMT2 = gdf$HTMT, se = se, lowerbound = lowerbound, upperbound = upperbound, time = tdelta, missing = NA)
  }else{
    print("ERROR")
  }
}

################################################################################
## jacknife: performs jacknife of a specified function (statisticfun) and 
# of a dataset (data), gives the confidence interval for the specified 
# siginificance level, mean and sd
################################################################################
jacknife <- function(data, statisticfun, ...,  alpha = 0.05)
{
  ind <- seq.int(from = 1, to = nrow(data))
  jack <- sapply(ind, function(x) statisticfun(data = data[-x, , drop = FALSE],...))  
  valid_jack <- jack[!is.na(jack)]
  lowerboundb <- unname(quantile(valid_jack, probs = alpha/2))
  upperboundb <- unname(quantile(valid_jack, probs = 1 - (alpha/2)))
  jackmean <- mean(valid_jack)
  jacksd <- sd(valid_jack)
  statcentered <- jackmean - valid_jack
  statcenteredsq <- (statcentered)^2
  statcenteredthree <- (statcentered)^3
  sumthree <- sum(statcenteredthree)
  sumsq <- sum(statcenteredsq)
  accelerator <- sumthree / (6*sumsq^(3/2))
  return(list(jacktab = jack, 
              mean = jackmean, 
              sd = jacksd, 
              accelerator = accelerator, 
              lowerbound = lowerboundb, 
              upperbound = upperboundb,
              missing = sum(is.na(jack))
  ))
}
################################################################################
## bootstrap: Performs bootstrap for the function (statisticfun) and of the dataset data.
## returns lowerbound and upperbound the mean and standard deviation of the statistic
################################################################################
bootstrap <- function(data, statisticfun, ...,  alpha = 0.05, nboot)
{
  boot <- sapply(1:nboot, function(x) statisticfun(data = dplyr::sample_n(data, nrow(data), replace = TRUE), ...))
  valid_boot <- boot[!is.na(boot)]
  lowerbound <- unname(quantile(valid_boot, probs = alpha/2))
  upperbound <- unname(quantile(valid_boot, probs = 1 - (alpha/2)))
  bootmean <- mean(valid_boot)
  bootsd <- sd(valid_boot)
  return(list(
    boot = valid_boot,
    lowerbound = lowerbound, 
    upperbound = upperbound,
    mean = bootmean, 
    se = bootsd,
    missing = sum(is.na(boot)),
    alpha = alpha
  ))
}
################################################################################
## bootbca: calculates the bootstrap bias corrected accelerated confidence intervals 
# with the help of the boot and jacknife functions. 
################################################################################
bootbca <- function(data, nboot, alpha = 0.05, statisticfun, ...){
  statistic <- statisticfun(data, ...)
  starttime <- Sys.time()
  boot <- bootstrap(data, nboot = nboot, alpha = alpha, statisticfun, ...)
  endtime <- Sys.time()
  tdeltaboot <- endtime - starttime
  z0 <- qnorm(p = mean(boot$boot < statistic), mean = 0, sd = 1)
  starttimebca <- Sys.time()
  jacknife <- jacknife(data, alpha = alpha, statisticfun, ...)
  acc <- jacknife$accelerator
  zalpha <- qnorm(p = alpha/2)
  zalpha2 <- qnorm(p = 1-(alpha/2))
  a1 <- pnorm(q = (z0 + (z0 + zalpha)/(1-acc*(z0+zalpha))), sd = 1, mean = 0)
  a2 <- pnorm(q = (z0 + (z0 + zalpha2)/(1-acc*(z0+zalpha2))), sd = 1, mean = 0)
  lowerbound <- unname(quantile(boot$boot, probs = a1))
  upperbound <- unname(quantile(boot$boot, probs = a2))
  endtimebca <- Sys.time()
  tdeltabootdca <- (endtimebca - starttimebca) + tdeltaboot
  return(list(boot = list(se = boot$se, 
                          lowerbound = boot$lowerbound, 
                          upperbound = boot$upperbound, 
                          time = tdeltaboot,
                          missing = boot$missing), 
              bootbca = list(se = boot$se, 
                             jack = jacknife$jacktab,
                             lowerbound = lowerbound, 
                             upperbound = upperbound, 
                             time = tdeltabootdca,
                             missing = jacknife$missing)))
} 
################################################################################
## wrapper: used in the simulation script. 
# uses the bootbca and deltamethod functions. 
# for each method the simulation scripts expects a lowerbound, upperbound and time. 
################################################################################
wrapper <- function(data, model, alpha, latent1, latent2, scale = scale, htmt2 = htmt2, nboot)
{
  delta <- deltamethod(data = data, model = model, alpha = alpha, 
                        latent1 = latent1, latent2 = latent2, 
                        scale = scale, htmt2 = htmt2)
  out <- bootbca(data = data, nboot = nboot, alpha = alpha, 
                 statisticfun = calchtmt, model = model_est, latent1 = latent1, 
                 latent2 <- latent2, scale = scale, htmt2 = htmt2)
  list(delta = delta, 
       boot = list(se = out$boot$se,
                   lowerbound = out$boot$lowerbound, 
                   upperbound = out$boot$upperbound,
                   missing = out$boot$missing,
                   time = out$boot$time), 
       bcaboot = list(se = out$bootbca$se,
                      lowerbound = out$bootbca$lowerbound, 
                      upperbound = out$bootbca$upperbound, 
                      time = out$bootbca$time,
                      missing = out$bootbca$missing
                      ))
}
################################################################################
## Defining the models
coefs <- 1
corr <- c(0.7, 0.8, 0.9, 0.95, 1)
param <- expand.grid(loading1 = coefs, loading2 = coefs, correlation = corr)
simModels <- foreach(i = 1:nrow(param), .combine = "rbind") %do%
  {
    simCommonFactor <- 
      paste(
        paste("xi_1 =~ ",param$loading1[i], "*x11 + ",param$loading1[i],"*x12 + ",param$loading1[i], "*x13"),"\n"
        , paste("xi_2 =~ ",param$loading2[i], "*x21 + ",param$loading2[i],"*x22 + ",param$loading2[i], "*x23"), "\n"
        , paste("xi_1 ~~ 1*xi_1 + ", param$correlation[i], "*xi_2"),"\n"
        , "xi_2 ~~ 1*xi_2 \n"
        , paste("x11 ~~", 0.6, "*x11 + 0*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x12 ~~", 0.5, "*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x13 ~~", 0.2, "*x13 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x21 ~~", 0.6, "*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x22 ~~", 0.5, "*x22 + 0*x23"),"\n"
        , paste("x23 ~~", 0.2, "*x23"), "\n"
        , paste("x11 ~ 0*1"), "\n"
        , paste("x12 ~ 0*1"), "\n"
        , paste("x13 ~ 0*1"), "\n"
        , paste("x21 ~ 0*1"), "\n"
        , paste("x22 ~ 0*1"), "\n"
        , paste("x23 ~ 0*1"), "\n"
      )
    save <- data.frame(
      loading_1 = param$loading1[i],
      loading_2 =  param$loading2[i],
      correlation = param$correlation[i],
      model = simCommonFactor
    )
    save
    #rm(save, simCommonFactor, i)
  }
simModels

rm(coefs, corr, param, simCommonFactor, save, i)
################################################################################
model_est<- '
              #  latent variables
                xi_1 =~ x11 + x12 + x13
                xi_2 =~ x21 + x22 + x23 
                
                xi_1 ~~ xi_2
              ' 