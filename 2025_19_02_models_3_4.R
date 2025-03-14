# nachbessern

coefs <- c(0.5, 0.7, 0.9)
loadings <- combn(coefs, 2)
#expand.grid function lieber als combn verwenden
expand.grid(loading_1 = coefs, loading_2 = coefs, correlation = c(0.5, 0.8, 0.9, 1))
correlation <- c(0.5, 0.8, 0.9, 1)
errorvar <- (1-loadings^2)
data_frame <- as.data.frame(matrix(nrow = , ncol = 10))
simModels <- foreach(i = 1:ncol(loadings), .combine = "rbind") %:% 
            foreach(j = c(0.5, 0.7, 0.8, 0.9, 1), .combine = "rbind") %do% 
            {
            simCommonFactor <- 
                           paste(
                             paste("xi_1 =~ ",loadings[1,i], "*x11 + ",loadings[1,i],"*x12 + ",loadings[1,i], "*x13 +",  loadings[1,i], "*x14"),"\n"
                           , paste("xi_2 =~ ",loadings[2,i], "*x21 + ",loadings[2,i],"*x22 + ",loadings[2,i], "*x23"), "\n"
                           , paste("xi_1 ~~ 1*xi_1 + ", j, "*xi_2"),"\n"
                           , "xi_2 ~~ 1*xi_2 \n"
                           , paste("x11 ~~", errorvar[1,i], "*x11 + 0*x12 + 0*x13 + 0*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
                           , paste("x12 ~~", errorvar[1,i], "*x12 + 0*x13 + 0*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
                           , paste("x13 ~~", errorvar[1,i], "*x13 + 0*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
                           , paste("x14 ~~", errorvar[1,i], "*x14 + 0*x21 + 0*x22 + 0*x23"),"\n"
                           , paste("x21 ~~", errorvar[2,i], "*x21 + 0*x22 + 0*x23"),"\n"
                           , paste("x22 ~~", errorvar[2,i], "*x22 + 0*x23"),"\n"
                           , paste("x23 ~~", errorvar[2,i], "*x23"), "\n"
                           )
           save <- data.frame(
                              loading_1 = loadings[1,i],
                              loading_2 = loadings[2,i],
                              correlation = j, 
                              model = simCommonFactor
                              )
           save
            }

