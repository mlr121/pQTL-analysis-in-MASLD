impute_data <- function(df, width = 0.3, downshift = 1.8) {
  # Ensure input is a data.frame
  df <- as.data.frame(df)
  
  # Only keep numeric columns for imputation
  numeric_cols <- names(df)[sapply(df, is.numeric)]
  
  # Make sure there are numeric columns
  if (length(numeric_cols) == 0) {
    stop("No numeric columns found in df.")
  }
  
  # Create logical imputation flags
  # impute_flags <- lapply(df[numeric_cols], function(x) !is.finite(x))
  #names(impute_flags) <- paste0(numeric_cols, "_impute")
  
  # Attach flags to the dataframe
  #df <- cbind(df, impute_flags)
  
  # Perform imputation
  set.seed(1)
  for (col in numeric_cols) {
    temp <- df[[col]]
    temp[!is.finite(temp)] <- NA
    sd_temp <- sd(temp, na.rm = TRUE)
    mean_temp <- mean(temp, na.rm = TRUE)
    
    # Define imputation parameters
    imp_sd <- width * sd_temp
    imp_mean <- mean_temp - downshift * sd_temp
    
    n_missing <- sum(is.na(temp))
    if (n_missing > 0) {
      temp[is.na(temp)] <- rnorm(n_missing, mean = imp_mean, sd = imp_sd)
    }
    df[[col]] <- temp
  }
  
  return(df)
}