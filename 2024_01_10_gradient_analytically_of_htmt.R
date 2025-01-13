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
  
  combos <- lapply(indicators, combn(,2))
  
  # Convert to data frame
  combos_df <- data.frame(
    var1 = combos[1,],
    var2 = combos[2,]
  )
  
}