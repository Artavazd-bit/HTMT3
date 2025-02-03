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


calc_grad_htmt_ana <- function(data, model, latent1, latent2, scale = TRUE){
  indicators <- extract_indicators(lv1 = latent1, lv2 = latent2, model_syntax = model)
  
  all_indicators <- unlist(indicators)
  
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
  
  cor_values$type[cor_values$col %in% unlist(indicators[1]) & cor_values$row %in% unlist(indicators[1])] <- "mono1"
  cor_values$type[cor_values$col %in% unlist(indicators[2]) & cor_values$row %in% unlist(indicators[2])] <- "mono2"
  cor_values$type[cor_values$col %in% unlist(indicators[1]) & cor_values$row %in% unlist(indicators[2])] <- "het"
  
  K_i <- length(unlist(indicators[1]))
  K_j <- length(unlist(indicators[2]))
  
  A = 1/(K_i*K_j) * sum(cor_values$val[cor_values$type == "het"])
  B = 2/(K_i*(K_i-1)) *  sum(cor_values$val[cor_values$type == "mono1"]) 
  C = 2/(K_j*(K_j-1)) *  sum(cor_values$val[cor_values$type == "mono2"]) 
  
  HTMT <- A / ((B*C)^(1/2))
  
  
  cor_values$gradient[cor_values$type == "het"] <- (1/(K_i*K_j) )/((B*C)^(1/2))
  cor_values$gradient[cor_values$type == "mono1"] <- A * 1/2 * 2/(K_i*(K_i-1)) * C * (B*C)^(-3/2) * -1
  cor_values$gradient[cor_values$type == "mono2"] <- A * 1/2 * 2/(K_j*(K_j-1)) * B * (B*C)^(-3/2) * -1
  
  list(output = cor_values, HTMT = HTMT)
} 



calc_grad_htmt2_ana <- function(data, model, latent1, latent2, scale = TRUE){
  indicators <- extract_indicators(lv1 = latent1, lv2 = latent2, model_syntax = model)
  
  all_indicators <- unlist(indicators)
  
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
  
  cor_values$type[cor_values$col %in% unlist(indicators[1]) & cor_values$row %in% unlist(indicators[1])] <- "mono1"
  cor_values$type[cor_values$col %in% unlist(indicators[2]) & cor_values$row %in% unlist(indicators[2])] <- "mono2"
  cor_values$type[cor_values$col %in% unlist(indicators[1]) & cor_values$row %in% unlist(indicators[2])] <- "het"
  
  K_i <- length(unlist(indicators[1]))
  K_j <- length(unlist(indicators[2]))
  
  A =  prod(cor_values$val[cor_values$type == "het"])^(1/(K_i*K_j))
  B =  prod(cor_values$val[cor_values$type == "mono1"])^(2/(K_i*(K_i-1))) 
  C =  prod(cor_values$val[cor_values$type == "mono2"])^(2/(K_j*(K_j-1))) 
  HTMT2 <- A / ((B*C)^(1/2))
  
  cor_values$gradient[cor_values$type == "het"] <- (1/(K_i*K_j)) * prod(cor_values$val[cor_values$type == "het"])^((1/(K_i*K_j))-1) * prod(cor_values$val[cor_values$type == "het"])/cor_values$val[cor_values$type == "het"] * 1/(sqrt((B*C)))
  cor_values$gradient[cor_values$type == "mono1"] <- A * 1/2 * (2/(K_i*(K_i-1))) * prod(cor_values$val[cor_values$type == "mono1"])^((2/(K_i*(K_i-1)))-1) * prod(cor_values$val[cor_values$type == "mono1"])/cor_values$val[cor_values$type == "mono1"] * C * (B*C)^(-3/2) * -1
  cor_values$gradient[cor_values$type == "mono2"] <- A * 1/2 * (2/(K_j*(K_j-1))) * prod(cor_values$val[cor_values$type == "mono2"])^((2/(K_j*(K_j-1)))-1) * prod(cor_values$val[cor_values$type == "mono2"])/cor_values$val[cor_values$type == "mono2"] * B * (B*C)^(-3/2) * -1
  
  list(output = cor_values, HTMT2 = HTMT2)
} 


calc_htmt <- function(data, model, latent1, latent2, scale = TRUE, htmt2 = FALSE){
  indicators <- extract_indicators(lv1 = latent1, lv2 = latent2, model_syntax = model)
  
  all_indicators <- unlist(indicators)
  
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
  
  cor_values$type[cor_values$col %in% unlist(indicators[1]) & cor_values$row %in% unlist(indicators[1])] <- "mono1"
  cor_values$type[cor_values$col %in% unlist(indicators[2]) & cor_values$row %in% unlist(indicators[2])] <- "mono2"
  cor_values$type[cor_values$col %in% unlist(indicators[1]) & cor_values$row %in% unlist(indicators[2])] <- "het"
  
  K_i <- length(unlist(indicators[1]))
  K_j <- length(unlist(indicators[2]))
  
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


