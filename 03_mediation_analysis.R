################################################################################
# Script: 03_mediation_analysis.R
# Description: Chr 2 pQTL mediation analysis for steatosis and protein targets.
################################################################################

################################################################################
# 1. Reproducibility & Environment Setup
################################################################################
set.seed(123)
options(stringsAsFactors = FALSE)

library(dplyr)
library(stringr)
library(here)
library(mediation)

# Source mediation helper function
source(here("functions", "mediation_helpers.R"))

# Simulation settings
sims_count <- 1000

################################################################################
# 2. Data Ingestion & Preprocessing
################################################################################
genotypes <- read.csv(here("data", "raw", "mediation_analysis_marker.csv")) %>%
  mutate(ID = str_replace_all(ID, "-", "_"))

data <- read.csv(here("data", "raw", "fixcol_ProMatrix_metadata.csv"))

steatosis <- read.csv(here("data", "raw", "pqtl_LIV_Stea.csv")) %>%
  mutate(ID = str_replace_all(ID, "-", "_"))

# Filter to Chr 2 pQTLs (Assumes pqtl_topLOD exists in environment or loaded)
chr2 <- pqtl_topLOD %>% filter(Chr == 2)

# Prepare protein matrix for Chr 2
data_chr2 <- data %>% filter(Gene.Symbol %in% chr2$Phenotype)
protein_matrix <- data_chr2
row.names(protein_matrix) <- protein_matrix$Gene.Symbol
protein_matrix$Gene.Symbol <- NULL
protein_matrix <- t(protein_matrix)

# Impute data (using custom/project impute_data function)
imputed_proteins <- impute_data(protein_matrix)
imputed_proteins$ID <- rownames(imputed_proteins)

# Merge datasets
med_data <- imputed_proteins %>%
  left_join(steatosis, by = "ID") %>%
  rename(steatosis = sum_all_vacuoles_percentage_24) %>%
  left_join(genotypes, by = "ID") %>%
  filter(!is.na(sex)) %>%
  rename_with(~str_remove_all(., "-"))

################################################################################
# 3. Identify Columns & Covariates
################################################################################
geno_cols <- grep("^Chr_2", colnames(med_data), value = TRUE)
protein_cols <- colnames(med_data)[!colnames(med_data) %in% c("ID", "steatosis", "sex", "Liver_g_24", "Liver_g_RelBWSac_24", geno_cols)]

covariate_cols <- c("Sex")

message("Found ", length(protein_cols), " protein targets and ", length(geno_cols), " genotype markers.")

################################################################################
# 4. Single Marker Mediation (Chr_2_125088490 -> Proteins -> Steatosis)
################################################################################
geno_fixed <- "Chr_2_125088490"

results_list_125 <- lapply(protein_cols, function(p) {
  run_mediation(
    df = med_data,
    treat_col = geno_fixed,
    mediator_col = p,
    outcome_col = "steatosis",
    covariate_cols = covariate_cols,
    sims = sims_count
  )
})

results_df_125 <- bind_rows(results_list_125)

# Filter to cis pQTLs on Chr 2 peak region (124Mb - 126Mb)
chr2_cis <- chr2 %>%
  filter(qtl_type == "cis" & Pos_bp > 124000000 & Pos_bp < 126000000)

cis_125 <- results_df_125 %>%
  filter(mediator %in% chr2_cis$Phenotype)

################################################################################
# 5. All Protein x Marker Mediation (Markers -> Proteins -> Steatosis)
################################################################################
results_list_all <- list()

for (g in geno_cols) {
  for (p in protein_cols) {
    res <- run_mediation(
      df = med_data,
      treat_col = g,
      mediator_col = p,
      outcome_col = "steatosis",
      covariate_cols = covariate_cols,
      sims = sims_count
    )
    if (!is.null(res)) {
      results_list_all[[paste(p, g, sep = "_")]] <- res
    }
  }
}

results_df <- bind_rows(results_list_all)

################################################################################
# 6. Fixed Mediator Analysis (Markers -> Ubr1 -> Target Proteins)
################################################################################
mediator_fixed <- "Ubr__1"
outcome_proteins <- setdiff(protein_cols, mediator_fixed)

results_ubr1 <- list()

for (g in geno_cols) {
  for (p in outcome_proteins) {
    res <- run_mediation(
      df = med_data,
      treat_col = g,
      mediator_col = mediator_fixed,
      outcome_col = p,
      covariate_cols = covariate_cols,
      sims = sims_count
    )
    if (!is.null(res)) {
      results_ubr1[[paste(mediator_fixed, p, g, sep = "_")]] <- res
    }
  }
}

results_ubr1_trans_df <- bind_rows(results_ubr1)

# Filter for significant indirect effects (full/partial mediation)
ubr1_sig_trans <- results_ubr1_trans_df %>%
  filter(p_ADE > 0.05 & p_ACME < 0.05)