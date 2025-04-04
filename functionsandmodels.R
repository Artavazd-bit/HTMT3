########################## Script for functions ################################

###################### Gradient Calculation HTMT and HTMT2 #####################
derivhtmt <- function(data, model, latent1, latent2, scale, htmt2){
  
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
  
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row %in% unlist(listind1)] <- "mono1"
  cor_values$type[cor_values$col %in% unlist(listind2) & cor_values$row %in% unlist(listind2)] <- "mono2"
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row %in% unlist(listind2)] <- "het"
  
  K_i <- length(unlist(listind1))
  K_j <- length(unlist(listind2))
  
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
    
    cor_values$gradient[cor_values$type == "het"] <- (1/(K_i*K_j)) * prod(cor_values$val[cor_values$type == "het"])^((1/(K_i*K_j))-1) * prod(cor_values$val[cor_values$type == "het"])/cor_values$val[cor_values$type == "het"] * 1/(sqrt((B*C)))
    cor_values$gradient[cor_values$type == "mono1"] <- A * 1/2 * (2/(K_i*(K_i-1))) * prod(cor_values$val[cor_values$type == "mono1"])^((2/(K_i*(K_i-1)))-1) * prod(cor_values$val[cor_values$type == "mono1"])/cor_values$val[cor_values$type == "mono1"] * C * (B*C)^(-3/2) * -1
    cor_values$gradient[cor_values$type == "mono2"] <- A * 1/2 * (2/(K_j*(K_j-1))) * prod(cor_values$val[cor_values$type == "mono2"])^((2/(K_j*(K_j-1)))-1) * prod(cor_values$val[cor_values$type == "mono2"])/cor_values$val[cor_values$type == "mono2"] * B * (B*C)^(-3/2) * -1
  }else{
    print("ERROR")
  }
  
list(output = cor_values, HTMT = HTMT)
} 

########################## Calculation HTMT/HTMT2 ##############################
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
  
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row %in% unlist(listind1)] <- "mono1"
  cor_values$type[cor_values$col %in% unlist(listind2) & cor_values$row %in% unlist(listind2)] <- "mono2"
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row %in% unlist(listind2)] <- "het"
  
  K_i <- length(unlist(listind1))
  K_j <- length(unlist(listind2))
  
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
}


############################# Calculation of Omega ############################# 
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
      products[, i, j] <- data_centered[, i] * data_centered[, j] # hier ändere ich die Objektklasse, drop = FALSE
      products[, j, i] <- products[, i, j]
    }
  }
  
  # Calculate covariances using vectorized operations
  for(idx1 in 1:nrow(indices)) {
    x <- indices[idx1, "row"] # nochmal indices anschauen
    y <- indices[idx1, "col"]
    for(idx2 in idx1:nrow(indices)) {
      z <- indices[idx2, "row"] # nochmal indices anschauen!!
      t <- indices[idx2, "col"]
      # Calculate fourth-order moments using pre-computed products
      omega_xyzt <- mean(products[, x, y] * products[, z, t]) - ( mean(products[, x, y]) * mean(products[, z, t]) )
      vc_r[idx1, idx2] <- omega_xyzt
      vc_r[idx2, idx1] <- omega_xyzt  # Matrix is symmetric
    }
  }
  return(vc_r)
}


# das sind einfach meine Daten
lavnacov <- function(data){
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
############################ Omega Correlation ################################
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



########################### Ultimate Function v3################################
doeverything <- function(data, model, alpha, latent1, latent2, scale, htmt2)
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
  zvalue <- (gdf$HTMT - 1)/se
  endtime <- Sys.time()
  tdelta <- endtime - starttime
  # here i want to test whether im in the lowest alpha percent cases. 
  ztest <- zvalue <  qnorm(p = alpha, mean = 0, sd = 1)
  
  if (htmt2 == FALSE){
    list(HTMT = gdf$HTMT, zval = zvalue, se = se, dec = ztest, time = tdelta)
  }else if (htmt2 == TRUE){
    list(HTMT2 = gdf$HTMT, zval = zvalue, se = se, dec = ztest, time = tdelta)
  }else{
    print("ERROR")
  }
}


simFun <- function(data, model, alpha, latent1, latent2, scale = scale, htmt2 = htmt2, seed, bootruns)
{
  delta <- doeverything(data = data, model = model, alpha = alpha, latent1 = latent1, latent2 = latent2, scale = scale, htmt2 = htmt2)
  
  
  set.seed(seed)
  starttime <- Sys.time()
  bootstrap <- boot(data_cfa, function(data, indices){calchtmt(data = data[indices,], model = model, latent1 = latent1, latent2 = latent2, scale = scale, htmt2 = htmt2)}, R = bootruns)
  bootstrap_htmt_1_se <- sd(bootstrap$t, na.rm = TRUE)
  #here its the other way, i want to test whether the upper bound is below zero
  bootupperbound <- quantile(bootstrap$t, na.rm = TRUE, probs = 1-alpha)
  endtime <- Sys.time()
  tdeltaboot <- endtime - starttime
  
  bootdec <- bootupperbound < 1
  
  list(delta = delta, boot = list(upperbound = bootupperbound, dec = bootdec, time = tdeltaboot))
}

############################# Models ###########################################


coefs <- c(0.5, 0.7, 0.9)
corr <- c(0.5, 0.8, 0.9, 1)
grid <- expand.grid(loading1 = coefs, loading2 = coefs, correlation = corr)
grid$errorvar1 <- (1-grid$loading1^2)
grid$errorvar2 <- (1-grid$loading2^2)
simModels_parallel <- foreach(i = 1:nrow(grid), .combine = "rbind") %do%
  {
    simCommonFactor <- 
      paste(
        paste("xi_1 =~ ",grid$loading1[i], "*x11 + ",grid$loading1[i],"*x12 + ",grid$loading1[i], "*x13 +",  grid$loading1[i], "*x14"),"\n"
        , paste("xi_2 =~ ",grid$loading2[i], "*x21 + ",grid$loading2[i],"*x22 + ",grid$loading2[i], "*x23"), "\n"
        , paste("xi_1 ~~ 1*xi_1 + ", grid$correlation[i], "*xi_2"),"\n"
        , "xi_2 ~~ 1*xi_2 \n"
        , paste("x11 ~~", grid$errorvar1[i], "*x11 + 0*x12 + 0*x13 + 0*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x12 ~~", grid$errorvar1[i], "*x12 + 0*x13 + 0*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x13 ~~", grid$errorvar1[i], "*x13 + 0*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x14 ~~", grid$errorvar1[i], "*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x21 ~~", grid$errorvar2[i], "*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x22 ~~", grid$errorvar2[i], "*x22 + 0*x23"),"\n"
        , paste("x23 ~~", grid$errorvar2[i], "*x23"), "\n"
      )
    save <- data.frame(
      loading_1 = grid$loading1[i],
      loading_2 =  grid$loading2[i],
      correlation = grid$correlation[i], 
      model = simCommonFactor
    )
    save
  }
simModels_parallel

################################################################################

coefs <- c(0.5, 0.7, 0.9)
corr <- c(0.5, 0.8, 0.9, 1)
grid <- expand.grid(loading1 = coefs, loading2 = coefs, correlation = corr)
grid$errorvar1 <- (1-grid$loading1^2)
grid$errorvar2 <- (1-grid$loading2^2)
simModels_tau <- foreach(i = 1:nrow(grid), .combine = "rbind") %do%
  {
    simCommonFactor <- 
      paste(
        paste("xi_1 =~ ",grid$loading1[i], "*x11 + ",grid$loading1[i],"*x12 + ",grid$loading1[i], "*x13 +",  grid$loading1[i], "*x14"),"\n"
        , paste("xi_2 =~ ",grid$loading2[i], "*x21 + ",grid$loading2[i],"*x22 + ",grid$loading2[i], "*x23"), "\n"
        , paste("xi_1 ~~ 1*xi_1 + ", grid$correlation[i], "*xi_2"),"\n"
        , "xi_2 ~~ 1*xi_2 \n"
        , paste("x11 ~~", 4, "*x11 + 0*x12 + 0*x13 + 0*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x12 ~~", 7, "*x12 + 0*x13 + 0*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x13 ~~", 5, "*x13 + 0*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x14 ~~", 3, "*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x21 ~~", 2, "*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x22 ~~", 1, "*x22 + 0*x23"),"\n"
        , paste("x23 ~~", 9, "*x23"), "\n"
      )
    save <- data.frame(
      loading_1 = grid$loading1[i],
      loading_2 =  grid$loading2[i],
      correlation = grid$correlation[i],
      model = simCommonFactor
    )
    save
  }
simModels_tau

################################################################################
coefs <- 1
corr <- c(0.7, 0.8, 0.9, 1)
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
      )
    save <- data.frame(
      loading_1 = param$loading1[i],
      loading_2 =  param$loading2[i],
      correlation = param$correlation[i],
      model = simCommonFactor
    )
    save
  }
simModels

################################################################################

model_est_old <- '
              #  latent variables
                xi_1 =~ x11 + x12 + x13 + x14
                xi_2 =~ x21 + x22 + x23 
                
                xi_1 ~~ xi_2
              ' 
###############################################################################

model_est<- '
              #  latent variables
                xi_1 =~ x11 + x12 + x13
                xi_2 =~ x21 + x22 + x23 
                
                xi_1 ~~ xi_2
              ' 


