extract_indicators <- function(lv1, lv2, model_syntax) {
  # Split the string into lines and remove empty lines
  lines <- strsplit(model_syntax, "\n")[[1]]
  lines <- trimws(lines[nchar(trimws(lines)) > 0])
  
  # Create a named list to store results
  result <- list()
  
  # Process each line
  for (line in lines) {
    # Remove comments and trim
    line <- trimws(gsub("#.*$", "", line))
    
    # Skip empty lines
    if (nchar(line) == 0) next
    
    # Split on =~
    parts <- strsplit(line, "=~")[[1]]
    
    # Get latent variable name (left side)
    lv_name <- trimws(parts[1])
    
    # If this is one of our target latent variables
    if (lv_name %in% c(lv1, lv2)) {
      # Get indicators (right side)
      indicators <- strsplit(trimws(parts[2]), "\\+")[[1]]
      indicators <- trimws(indicators)
      
      # Store in result list
      result[[lv_name]] <- indicators
    }
  }
  
  return(result)
}


calc_grad_htmt_anal <- function(data, model, latent1, latent2){
  indicators <- extract_indicators(lv1 = latent1, lv2 = latent2, model_syntax = model_est)
  
  all_indicators <- unlist(indicators)
  
  subset_data <- data[, all_indicators]
  
  cor_subset_data <- cor(subset_data)
  
  cor_subset_data_mono1 <- cor_subset_data[unlist(indicators[1], FALSE), unlist(indicators[1], FALSE)]
  cor_subset_data_mono2 <- cor_subset_data[unlist(indicators[2], FALSE), unlist(indicators[2], FALSE)]
  cor_subset_data_het <- cor_subset_data[unlist(indicators[2], FALSE), unlist(indicators[1], FALSE)]
  
  cor_subset_data_mono1_val <- cor_subset_data_mono1[lower.tri(cor_subset_data_mono1)]
  cor_subset_data_mono2_val <- cor_subset_data_mono2[lower.tri(cor_subset_data_mono2)]
  
  K_i <- length(unlist(indicators[1]))
  K_j <- length(unlist(indicators[2]))
  
  A = 1/(K_i*K_j) * sum(cor_subset_data_het)
  B = 2/(K_i*(K_i-1)) *  sum(cor_subset_data_mono1_val) 
  C = 2/(K_j*(K_j-1)) *  sum(cor_subset_data_mono2_val) 
  
  HTMT <- A / ((B*C)^(1/2))
  
  het_grad <- (1/(K_i*K_j) )/((B*C)^(1/2))
  
  mono1_grad <- A * 1/2 * 2/(K_i*(K_i-1)) * C * (B*C)^(-3/2) * -1
  mono2_grad <-  A * 1/2 * 2/(K_j*(K_j-1)) * B * (B*C)^(-3/2) * -1
  
  c(het_grad,
  mono1_grad,
  mono2_grad)
} 



calc_grad_htmt2_anal <- function(data, model, latent1, latent2){
  indicators <- extract_indicators(lv1 = latent1, lv2 = latent2, model_syntax = model_est)
  
  all_indicators <- unlist(indicators)
  
  subset_data <- data[, all_indicators]
  
  cor_subset_data <- cor(subset_data)
  
  cor_subset_data_mono1 <- cor_subset_data[unlist(indicators[1], FALSE), unlist(indicators[1], FALSE)]
  cor_subset_data_mono2 <- cor_subset_data[unlist(indicators[2], FALSE), unlist(indicators[2], FALSE)]
  cor_subset_data_het <- cor_subset_data[unlist(indicators[2], FALSE), unlist(indicators[1], FALSE)]
  
  cor_subset_data_mono1_val <- cor_subset_data_mono1[lower.tri(cor_subset_data_mono1)]
  cor_subset_data_mono2_val <- cor_subset_data_mono2[lower.tri(cor_subset_data_mono2)]
  
  K_i <- length(unlist(indicators[1]))
  K_j <- length(unlist(indicators[2]))
  
  A =  prod(cor_subset_data_het)^(1/(K_i*K_j))
  B =  prod(cor_subset_data_mono1_val)^(2/(K_i*(K_i-1))) 
  C =  prod(cor_subset_data_mono2_val)^(2/(K_j*(K_j-1))) 
  HTMT2 <- A / ((B*C)^(1/2))
  
  all_val <- c(mon1o = cor_subset_data_mono1_val, het = cor_subset_data_het, mon2o =cor_subset_data_mono2_val)
  cols <-  gsub("\\d+$", "", names(all_val))
  
  gradient_htmt2 <- rep(NA, length(all_val))
  #print(cols)
  i <- 1
  for(x in all_val){
    #print(x)
    if(cols[i] == "het"){
      #print(i)
      grad_A = (1/(K_i*K_j)) * prod(cor_subset_data_het)^((1/(K_i*K_j))-1) * prod(cor_subset_data_het)/x
      gradient_htmt2[i] <- grad_A * 1/(sqrt((B*C)))
      #print(gradient_htmt2[i])
    }else if(cols[i] == "mon1o"){
      #print(i)
      grad_B = (2/(K_i*(K_i-1))) * prod(cor_subset_data_mono1_val)^((2/(K_i*(K_i-1)))-1) * prod(cor_subset_data_mono1_val)/x
      gradient_htmt2[i] <- A * 1/2 * grad_B * C * (B*C)^(-3/2) * -1
      #print(gradient_htmt2[i])
    }else if(cols[i] == "mon2o"){
      #print(i)
      grad_C = (2/(K_j*(K_j-1))) * prod(cor_subset_data_mono2_val)^((2/(K_j*(K_j-1)))-1) * prod(cor_subset_data_mono2_val)/x
      gradient_htmt2[i] <- A * 1/2 * grad_C * B * (B*C)^(-3/2) * -1
      #print(gradient_htmt2[i])
    }
    else{
      #print(i)
      gradient_htmt2[i] <- NA
    }
    i <- i + 1
  }
  #lapply(1:length(cor_subset_data_het), function(i) prod(cor_subset_data_het[-i]))
  gradient_htmt2
} 

grad_htmtl2_anal <- calc_grad_htmt2_anal(data = data_cfa, model = model_est, latent1 = "xi_1", latent2 = "xi_2")



