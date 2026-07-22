#load data ----
genotypes <- read.csv("raw_data\\mediation_analysis_marker.csv", stringsAsFactors = F)

data <- read.csv(file="raw_data\\fixcol_ProMatrix_metadata.csv", header = TRUE, stringsAsFactors = FALSE)

steatosis <- read.csv("raw_data\\pqtl_LIV_Stea.csv", stringsAsFactors = F)
steatosis$ID <- gsub("-","_", steatosis$ID)
library(dplyr)

chr2 <- pqtl_topLOD %>%
  filter(Chr == 2)

data_chr2 <- data %>%
  filter(Gene.Symbol %in% chr2$Phenotype)

protein_matrix <- data_chr2
row.names(protein_matrix) <- protein_matrix$Gene.Symbol
protein_matrix$Gene.Symbol <-NULL
protein_matrix <- t(protein_matrix)

imputed_proteins <- impute_data(protein_matrix)
imputed_proteins$ID <- rownames(imputed_proteins)

imputed_proteins <- imputed_proteins %>%
  left_join(steatosis, by = c("ID" = "ID")) %>%
  rename(steatosis = sum_all_vacuoles_percentage_24)

genotypes$ID <- gsub("-","_", genotypes$ID)
med_data <- imputed_proteins %>%
  left_join(genotypes, by = c("ID" = "ID")) %>%
  filter(!is.na(sex)) %>%
  rename_with(~str_remove_all(., "-"))

# Identify column groups ----
geno_cols <- grep("^Chr_2", colnames(med_data), value = TRUE)
protein_cols <- colnames(med_data)[!colnames(med_data) %in% c("ID", "steatosis", "sex", "Liver_g_24", "Liver_g_RelBWSac_24", geno_cols)]

# Confirm
length(protein_cols)
length(geno_cols)
head(protein_cols)
head(geno_cols)

  #founder_levels <- c("B6-129", "CAST-129", "B6-PWK", "CAST-PWK")

covariate_cols <- c("Sex")

#single marker mediation ----
geno_fixed <- "Chr_2_125088490"

min_samples <- floor(0.5 * nrow(med_data))
results_list <- list()

for (p in protein_cols) {
  
  relevant_cols <- c(geno_fixed, p, "steatosis", covariate_cols)
  df <- med_data[complete.cases(med_data[, relevant_cols]), relevant_cols]
  
  # QC checks
  if (nrow(df) < min_samples) {
    message("❌ Skipping ", p, ": only ", nrow(df), " complete cases (<50%)")
    next
  }
  if (length(unique(df[[geno_fixed]])) < 2) next
  if (length(unique(df[[p]])) < 2) next
  
  message("✅ Running mediation for ", p, " × ", geno_fixed,
          " with ", nrow(df), " complete cases")
  
  # Build formulas — backtick-wrap all names to handle invalid R identifiers
  cov_str <- if (length(covariate_cols) > 0)
    paste("+", paste(paste0("`", covariate_cols, "`"), collapse = " + "))
  else ""
  
  med_formula <- as.formula(paste0("`", p, "` ~ `", geno_fixed, "`", cov_str))
  out_formula <- as.formula(paste0("steatosis ~ `", p, "` + `", geno_fixed, "`", cov_str))
  
  # Fit models
  med_model <- lm(med_formula, data = df)
  out_model <- lm(out_formula, data = df)
  
  # Mediation
  med_out <- tryCatch({
    mediate(
      model.m  = med_model,
      model.y  = out_model,
      treat    = geno_fixed,
      mediator = p,
      boot     = TRUE,
      sims     = sims
    )
  }, error = function(e) {
    message("⚠️ Error for ", p, ": ", e$message)
    NULL
  })
  
  # Save results
  if (!is.null(med_out)) {
    results_list[[p]] <- data.frame(
      mediator      = p,
      genotype      = geno_fixed,
      ACME          = med_out$d0,
      ADE           = med_out$z0,
      total_effect  = med_out$tau.coef,
      prop_mediated = med_out$n0,
      p_ACME        = med_out$d0.p,
      p_ADE         = med_out$z0.p
    )
  }
}

results_df_125 <- bind_rows(results_list)

##filter to cis on chr 2 peak ----
chr2_cis <- chr2 %>%
  filter(qtl_type == "cis") %>%
  filter(pos.bp.> 124000000 & pos.bp.< 126000000)
cis_125 <- results_df_125 %>%
  filter(mediator %in% chr2_cis$Phenotype)


# all protein x gene mediation ----

# Set a minimum number of complete cases (or max missingness)
min_samples <- floor(0.5 * nrow(med_data))  # keep pairs with >=50% non-missing

results_list <- list()

for (g in geno_cols) {
  for (p in protein_cols) {
    # Subset to complete cases
    relevant_cols <- c(g, p, "steatosis", covariates)
    df <- med_data[complete.cases(med_data[, relevant_cols]), relevant_cols]
    
    # Skip if too few complete cases
    if (nrow(df) < min_samples) {
      message("❌ Skipping ", p, " × ", g, ": only ", nrow(df), " complete cases (<50%)")
      next
    }
    
    # Skip if genotype has <2 levels
    if (length(unique(df[[g]])) < 2) {
      message("❌ Skipping ", p, " × ", g, ": genotype has <2 levels")
      next
    }
    
    # Skip if mediator has <2 levels
    if (length(unique(df[[p]])) < 2) {
      message("❌ Skipping ", p, " × ", g, ": mediator has <2 levels")
      next
    }
    
    message("✅ Running mediation for ", p, " × ", g, " with ", nrow(df), " complete cases")
    
    # Build formulas
    cov_str <- if (length(covariates) > 0) paste("+", paste(covariates, collapse = " + ")) else ""
    med_formula <- as.formula(paste(p, "~", g, cov_str))
    out_formula <- as.formula(paste("steatosis ~", p, "+", g, cov_str))
    
    # Fit models
    med_model <- lm(med_formula, data = df)
    out_model <- lm(out_formula, data = df)
    
    # Run mediation
    med_out <- tryCatch({
      mediate(
        model.m = med_model,
        model.y = out_model,
        treat = g,
        mediator = p,
        boot = TRUE,
        sims = sims
      )
    }, error = function(e) {
      message("⚠️ Error for ", p, " × ", g, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(med_out)) {
      results_list[[paste(p, g, sep = "_")]] <- data.frame(
        genotype = g,
        mediator = p,
        ACME = med_out$d0,
        ADE = med_out$z0,
        total_effect = med_out$tau.coef,
        prop_mediated = med_out$n0,
        p_ACME = med_out$d0.p,
        p_ADE = med_out$z0.p
      )
    }
  }
}

results_df <- bind_rows(results_list)

# Mediator fixed ----
mediator_fixed <- "Ubr__1"

# Outcomes: all proteins except Ubr1
outcome_proteins <- setdiff(protein_cols, mediator_fixed)

results_ubr1 <- list()
min_samples <- floor(0.5 * nrow(med_data))  # 50% complete cases threshold

for (g in geno_cols) {
  for (p in outcome_proteins) {
    relevant_cols <- c(g, mediator_fixed, p, covariate_cols)
    df <- med_data[complete.cases(med_data[, relevant_cols]), relevant_cols]
    
    # Skip if too few complete cases
    if (nrow(df) < min_samples) {
      message("❌ Skipping ", mediator_fixed, " → ", p, " × ", g, ": only ", nrow(df), " complete cases (<50%)")
      next
    }
    
    # Skip if genotype has <2 levels
    if (length(unique(df[[g]])) < 2) next
    
    # Skip if mediator or outcome has <2 levels
    if (length(unique(df[[mediator_fixed]])) < 2) next
    if (length(unique(df[[p]])) < 2) next
    
    message("✅ Running mediation for ", mediator_fixed, " → ", p, " × ", g, " with ", nrow(df), " complete cases")
    
    cov_str <- if (length(covariate_cols) > 0) paste("+", paste(covariate_cols, collapse = " + ")) else ""
    med_formula <- as.formula(paste(mediator_fixed, "~", g, cov_str))
    out_formula <- as.formula(paste(p, "~", mediator_fixed, "+", g, cov_str))
    
    # Fit models
    med_model <- lm(med_formula, data = df)
    out_model <- lm(out_formula, data = df)
    
    # Run mediation
    med_out <- tryCatch({
      mediate(
        model.m = med_model,
        model.y = out_model,
        treat = g,
        mediator = mediator_fixed,
        boot = TRUE,
        sims = sims
      )
    }, error = function(e) {
      message("⚠️ Error for ", mediator_fixed, " → ", p, " × ", g, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(med_out)) {
      results_ubr1[[paste(mediator_fixed, p, g, sep = "_")]] <- data.frame(
        genotype = g,
        mediator = mediator_fixed,
        outcome = p,
        ACME = med_out$d0,
        ADE = med_out$z0,
        total_effect = med_out$tau.coef,
        prop_mediated = med_out$n0,
        p_ACME = med_out$d0.p,
        p_ADE = med_out$z0.p
      )
    }
  }
}

results_ubr1_trans_df <- bind_rows(results_ubr1)

ubr1_sig_trans <- results_ubr1_trans_df %>%
  filter(p_ADE > 0.05 & p_ACME < 0.05)