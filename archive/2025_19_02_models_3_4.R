library(foreach)

###############################################################################################################

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

###############################################################################################################

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
        , paste("x11 ~~", 0.6, "*x11 + 0*x12 + 0*x13 + 0*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x12 ~~", 0.7, "*x12 + 0*x13 + 0*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x13 ~~", 0.8, "*x13 + 0*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x14 ~~", 0.9, "*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x21 ~~", 0.7, "*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x22 ~~", 0.9, "*x22 + 0*x23"),"\n"
        , paste("x23 ~~", 0.5, "*x23"), "\n"
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

###############################################################################################################



