################################################################################
# Script: mediation_helpers.R
# Description: Helper function to fit mediation models safely.
################################################################################

library(mediation)
library(dplyr)

run_mediation <- function(df, treat_col, mediator_col, outcome_col, covariate_cols = NULL, min_prop = 0.5, sims = 1000) {
  relevant_cols <- c(treat_col, mediator_col, outcome_col, covariate_cols)
  
  # Ensure columns exist
  if (!all(relevant_cols %in% colnames(df))) {
    missing <- relevant_cols[!relevant_cols %in% colnames(df)]
    message("⚠️ Missing columns in dataframe: ", paste(missing, collapse = ", "))
    return(NULL)
  }
  
  # Subset complete cases
  sub_df <- df[complete.cases(df[, relevant_cols]), relevant_cols, drop = FALSE]
  
  # QC checks
  if (nrow(sub_df) < floor(min_prop * nrow(df))) return(NULL)
  if (length(unique(sub_df[[treat_col]])) < 2) return(NULL)
  if (length(unique(sub_df[[mediator_col]])) < 2) return(NULL)
  if (length(unique(sub_df[[outcome_col]])) < 2) return(NULL)
  
  # Construct formulas safely with backticks
  cov_str <- if (length(covariate_cols) > 0) {
    paste("+", paste0("`", covariate_cols, "`", collapse = " + "))
  } else {
    ""
  }
  
  med_formula <- as.formula(paste0("`", mediator_col, "` ~ `", treat_col, "` ", cov_str))
  out_formula <- as.formula(paste0("`", outcome_col, "` ~ `", mediator_col, "` + `", treat_col, "` ", cov_str))
  
  # Fit linear models
  med_model <- lm(med_formula, data = sub_df)
  out_model <- lm(out_formula, data = sub_df)
  
  # Run causal mediation
  med_out <- tryCatch({
    mediation::mediate(
      model.m  = med_model,
      model.y  = out_model,
      treat    = treat_col,
      mediator = mediator_col,
      boot     = TRUE,
      sims     = sims
    )
  }, error = function(e) {
    message("⚠️ Error for ", treat_col, " -> ", mediator_col, " -> ", outcome_col, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(med_out)) return(NULL)
  
  # Return structured output
  data.frame(
    genotype      = treat_col,
    mediator      = mediator_col,
    outcome       = outcome_col,
    ACME          = med_out$d0,
    ADE           = med_out$z0,
    total_effect  = med_out$tau.coef,
    prop_mediated = med_out$n0,
    p_ACME        = med_out$d0.p,
    p_ADE         = med_out$z0.p,
    stringsAsFactors = FALSE
  )
}