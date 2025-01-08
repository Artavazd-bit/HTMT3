####################################################
# This function calculates the variance covariance matrix of a correlation matrix


calculate_corr_cov <- function(data){
  ####################################################################################
  ## Berechnung der Varianz-Covarianz-Matrix der Korrelationskoeffizienten gilt nur asymptotisch Dykstra(2013) S.11
  ## nach Dykstra(2013) Gleichung 27/28 und Isserlis(2019) Gleichung 21 + Anmerkungen von Flo
  ####################################################################################
  size_n <- ((ncol(data)*ncol(data)) - ncol(data))/2
  vc_r <- matrix(0 , nrow = size_n , ncol = size_n)
  cor_data <- cor(data)
  cnter_y <- 1
  cnter_t <- 1
  for(x in 1:(ncol(data)-1)){
    bar_x <- mean(data[,x])
    sigma_x <- sd(data[,x])
    for(y in (x+1):ncol(data)){
      bar_y <- mean(data[,y])
      sigma_y <- sd(data[,y])
      for(z in 1:(ncol(data)-1)){
        bar_z <- mean(data[,z])
        sigma_z <- sd(data[,z])
        for(t in (z+1):ncol(data)){
          bar_t <- mean(data[,t])
          sigma_t <- sd(data[,t])
          
          mu_xyzt_temp <- 0
          mu_xxzt_temp <- 0 
          mu_yyzt_temp <- 0
          mu_xyzz_temp <- 0 
          mu_xytt_temp <- 0
          
          mu_xxtt_temp <- 0
          mu_xxzz_temp <- 0
          
          mu_yytt_temp <- 0
          mu_yyzz_temp <- 0
          
          for(i in 1:nrow(data)){
            mu_xyzt_temp <- mu_xyzt_temp + (data[i,x] - bar_x) * (data[i,y] - bar_y) * (data[i,z] - bar_z) * (data[i,t] - bar_t)
            
            mu_xxzt_temp <- mu_xxzt_temp + (data[i,x] - bar_x) * (data[i,x] - bar_x) * (data[i,z] - bar_z) * (data[i,t] - bar_t)
            mu_yyzt_temp <- mu_yyzt_temp + (data[i,y] - bar_y) * (data[i,y] - bar_y) * (data[i,z] - bar_z) * (data[i,t] - bar_t)
            mu_xyzz_temp <- mu_xyzz_temp + (data[i,x] - bar_x) * (data[i,y] - bar_y) * (data[i,z] - bar_z) * (data[i,z] - bar_z)
            mu_xytt_temp <- mu_xytt_temp + (data[i,x] - bar_x) * (data[i,y] - bar_y) * (data[i,t] - bar_t) * (data[i,t] - bar_t)
            
            mu_xxtt_temp <- mu_xxtt_temp + (data[i,x] - bar_x) * (data[i,x] - bar_x) * (data[i,t] - bar_t) * (data[i,t] - bar_t)
            mu_xxzz_temp <- mu_xxzz_temp + (data[i,x] - bar_x) * (data[i,x] - bar_x) * (data[i,z] - bar_z) * (data[i,z] - bar_z)
            
            mu_yytt_temp <- mu_yytt_temp + (data[i,y] - bar_y) * (data[i,y] - bar_y) * (data[i,t] - bar_t) * (data[i,t] - bar_t)
            mu_yyzz_temp <- mu_yyzz_temp + (data[i,y] - bar_y) * (data[i,y] - bar_y) * (data[i,z] - bar_z) * (data[i,z] - bar_z)
          }
          
          mu_xyzt <- 1/nrow(data) * mu_xyzt_temp
          
          mu_xxzt <- 1/nrow(data) * mu_xxzt_temp
          
          mu_yyzt <- 1/nrow(data) * mu_yyzt_temp
          
          mu_xyzz <- 1/nrow(data) * mu_xyzz_temp
          
          mu_xytt <- 1/nrow(data) * mu_xytt_temp
          
          
          mu_xxtt <- 1/nrow(data) * mu_xxtt_temp
          mu_xxzz <- 1/nrow(data) * mu_xxzz_temp
          
          mu_yytt <- 1/nrow(data) * mu_yytt_temp
          mu_yyzz <- 1/nrow(data) * mu_yyzz_temp
          
          vc_r[cnter_y, cnter_t]<- (mu_xyzt - 1/2*cor_data[x,y]*(mu_xxzt + mu_yyzt) - 1/2*cor_data[z,t]*(mu_xyzz + mu_xytt) 
                                    + 1/4*cor_data[x,y]*cor_data[z,t]*(mu_xxzz+mu_xxtt+mu_yyzz+mu_yytt))/nrow(data)
          cnter_t <- cnter_t + 1
        }
      }
      cnter_y <- cnter_y + 1
      cnter_t <- 1
    }
  }
  return(vc_r)
}

library(stringr)

calc_gradient <- function(data, sim_runs, n, jj, model_est, latent1_index = 1, latent2_index = 2){
  
  latent_vars <- str_extract_all(model_est, "xi_[0-9]+")[[1]]
  x_vars <- str_extract_all(model_est, "x\\d+")[[1]]
  cnt_x_vars <- length(x_vars)
  cnt_cor_val <- (cnt_x_vars*cnt_x_vars - cnt_x_vars)/2
  
  cor_data <- cor(data)
  htmt2_base <- semTools::htmt(model = model_est,
                          data =  NULL, 
                          sample.cov = cor_data,
                          htmt2 = TRUE
  )
  htmt_base <-  semTools::htmt(model = model_est,
                          data =  NULL, 
                          sample.cov = cor_data,
                          htmt2 = FALSE
  )
  
  # für den Gradienten der pfadestimates
  gr = matrix(0 , nrow = size_n_ , ncol = 1)
  gdr = matrix(0 , nrow = size_n_ , ncol = 1)
  # Path estimates
  dlta_bt_htmt2 = matrix(0 , nrow = cnt_x_vars , ncol = cnt_cor_val)
  dlta_bt_htmt = matrix(0 , nrow = cnt_x_vars, ncol = cnt_cor_val)
  # für die h - Methode
  list_dh <- c(0.00001, 0.0001, 0.001, 0.01, 0.1)
  #Übersicht für verschiedene dh
  gradient_htmt <- data.frame(matrix(0, nrow = cnt_cor_val, ncol = length(list_dh)))
  colnames(gradient_htmt) <- list_dh
  
  gradient_htmt2 <-data.frame(matrix(0, nrow = cnt_cor_val, ncol = 1))
  colnames(gradient_htmt2) <- list_dh
  
  
  # Counter, sodass jeder Korrelationskoeffizient durchgegangen wird
  cnter = 1
  for(dh in list_dh){
    #print(dh)
    for (i in 1:(ncol(data)-1)) {
      for (j in (i+1):ncol(data)) {
        cor_data <- cor(data)
        cor_data[i,j] = cor_data[i,j] + dh
        cor_data[j,i] = cor_data[i,j]
        set.seed(50+jj+sim_runs+n)
        
        htmt2 <- semTools::htmt(model = model_est,
                       data =  NULL, 
                       sample.cov = cor_data,
                       htmt2 = TRUE
                        )
        htmt <-  semTools::htmt(model = model_est,
                                data =  NULL, 
                                sample.cov = cor_data,
                                htmt2 = FALSE
        )
        # in h-Methode: f(x+h):
        bt_htmt2 = htmt2[latent_vars[latent1_index], latent_vars[latent1_index]]
        bt_htmt = htmt[latent_vars[latent1_index], latent_vars[latent1_index]]
        # f(x+h) - f(x) / h
        dlta_bt_htmt2[, cnter] = (bt_htmt2 - htmt2_base)/dh
        dlta_bt_htmt[, cnter] = (bt_htmt - htmt_base)/dh
        
        cnter = cnter + 1
      }
    }
    cnter <- 1
    gradient_htmt[, paste0(dh)]  <- dlta_bt_htmt
    gradient_htmt2[, paste0(dh)]  <- dlta_bt_htmt2
    Gbt = matrix(0 , nrow = Anzahl_path_estimate , ncol = size_n_)
    Gbt2 = matrix(0 , nrow = Anzahl_path_estimate , ncol = size_n_)
  }
  
  Gbt_x1 <- t(as.matrix(Gbt_overview_delta[,"1e-05"]))
  Gbt_x2 <- t(as.matrix(Gbt_overview_delta[,"1e-05"]))
}


