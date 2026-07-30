################################################################################
# Script: 04_ptm_ubiquitination_analysis.R
# Description: Ubiquitination site normalization, N-end rule classification,
#              Ubr1 GAM associations, attenuation modeling, and GSEA.
################################################################################

################################################################################
# 1. Reproducibility & Environment Setup
################################################################################
set.seed(123)
options(stringsAsFactors = FALSE)

# Core Libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(plotly)
library(purrr)
library(broom)
library(here)

# Genomics & Functional Libraries
library(Biostrings)
library(msigdbr)
library(fgsea)

# Source existing helper functions
source(here("functions", "gam_helpers.R"))

################################################################################
# 2. Data Ingestion & Integration
################################################################################
# Load PTM (Ubiquitination) Data
ptm <- read.delim(here("data", "raw", "simple_UbReport.txt"), stringsAsFactors = FALSE)

ub_matrix <- ptm %>%
  mutate(
    PG.Genes = paste0(Gene, EG.PrecursorId),
    PEP.IsProteinGroupSpecific = NULL,
    EG.PrecursorId = NULL,
    PG.FastaHeaders = NULL,
    EG.ProteinPTMLocations = NULL,
    Gene = NULL,
    PEP.StrippedSequence = NULL
  )

rownames(ub_matrix) <- ub_matrix$PG.Genes
ub_matrix$PG.Genes <- NULL

# Transpose and format samples
ub <- as.data.frame(t(ub_matrix)) %>%
  rename_with(~ str_remove(., "^X")) %>%
  mutate(across(everything(), as.numeric))

# Keep mapping table for sequence lookup
ub_ss <- ptm %>%
  select(Gene, EG.PrecursorId, PEP.StrippedSequence) %>%
  mutate(PG.Genes = paste0(Gene, EG.PrecursorId))

# Load Proteomics & Phenotypes
proteomics_raw <- read.csv(here("data", "raw", "fixcol_ProMatrix_metadata.csv"))
pheno <- read.csv(here("data", "raw", "F2_phenotypes.csv"))

proteomics_raw <- proteomics_raw %>%
  mutate(Gene.Symbol = ifelse(is.na(Gene.Symbol) | Gene.Symbol == "", protein.groups, Gene.Symbol))

proteomics_clean <- proteomics_raw %>%
  select(-c(Gene.Symbol, protein.groups, PG.FastaHeaders))

proteomics <- as.data.frame(t(proteomics_clean))
rownames(proteomics) <- str_remove_all(rownames(proteomics), "_")
colnames(proteomics) <- make.unique(proteomics_raw$Gene.Symbol)
proteomics$ID <- rownames(proteomics)

pheno <- pheno %>% mutate(ID = str_remove_all(ID, "-"))
non_traits <- c("ID", "Strain", "Generation", "Sex")
pheno_cols <- colnames(pheno)[!(colnames(pheno) %in% non_traits)]
pheno_cols <- pheno_cols[sapply(pheno[pheno_cols], is.numeric)]

all_pheno <- pheno %>% left_join(proteomics, by = "ID")

# Stitch Ubiquitination sites with phenotypes & total proteomics
ub$ID <- rownames(ub)
ubsites <- ub %>% left_join(all_pheno, by = "ID")

################################################################################
# 3. Ubiquitination Site Normalization by Total Protein Abundance
################################################################################
# Log2 transform expression and site columns
ub_log <- ubsites %>%
  mutate(across(where(is.numeric) & !all_of(pheno_cols), log2))

expr_cols <- setdiff(colnames(ub_log), pheno_cols)
ub_cols <- expr_cols[str_detect(expr_cols, "_")]
protein_cols <- setdiff(expr_cols, ub_cols)

ub_genes <- sub("_.*$", "", ub_cols)
has_protein <- ub_genes %in% protein_cols

ub_cols_matched <- ub_cols[has_protein]
ub_genes_matched <- ub_genes[has_protein]

ub_norm <- ub_log

# Subtract total log2 protein abundance from log2 Ub site abundance
for (i in seq_along(ub_cols_matched)) {
  site_col <- ub_cols_matched[i]
  prot_col <- ub_genes_matched[i]
  ub_norm[[site_col]] <- ub_log[[site_col]] - ub_log[[prot_col]]
}

################################################################################
# 4. Exploratory Visualization & GAM Associations
################################################################################
# Protein level plot
p_protein <- ggplot(all_pheno, aes(x = Ubr1, y = Plin2, alpha = Steatosis_perc_24)) +
  geom_point(aes(color = Sex)) +
  geom_smooth(method = "gam") +
  theme_minimal()

# Site level plot
target_site <- "Mcmbp_LQHINPLLPTC.Carbamidomethyl..C..LNK.GlyGly..K..EESR_.4"
if (target_site %in% colnames(ub_norm)) {
  p_site <- ggplot(ub_norm, aes(x = Ubr1, y = .data[[target_site]], alpha = Steatosis_perc_24)) +
    geom_point(aes(color = Sex), size = 3) +
    theme_minimal()
  ggplotly(p_site)
}

# Run GAM & Kendall Tau calculations via helper functions
ubr1_results      <- compute_ubr1_associations(ub_norm, ubr1_col = "Ubr1")
stea_tau_results  <- compute_ubr1_associations(ub_norm, ubr1_col = "Steatosis_perc_24")

# Combine Ub site associations
ubnorm_tau <- stea_tau_results %>%
  filter(Protein %in% ub_cols) %>%
  mutate(Gene = sub("_.*", "", Protein)) %>%
  left_join(ubr1_results, by = "Protein", suffix = c("_stea", "_ubr1"))

ubnorm_tau_sig <- ubnorm_tau %>% filter(p_value_ubr1 <= 0.05)

p_ubtau <- ggplot(ubnorm_tau, aes(x = tau_stea, y = tau_ubr1)) +
  geom_point(aes(text = Protein)) +
  theme_minimal() +
  labs(x = "Tau (Steatosis)", y = "Tau (UBR1)")

ggplotly(p_ubtau)

################################################################################
# 5. N-Terminal Residue Extraction & N-End Rule Classification
################################################################################
fasta_files <- c(
  here("data", "raw", "Mus_musculus_c57bl6nj.C57BL_6NJ_v1.pep.all.fasta"),
  here("data", "raw", "Mus_musculus_casteij.CAST_EiJ_v1.pep.all.fasta"),
  here("data", "raw", "Mus_musculus_pwkphj.PWK_PhJ_v1.pep.all.fasta"),
  here("data", "raw", "Mus_musculus_129s1svimj.129S1_SvImJ_v1.pep.all.fasta")
)

# Filter existing files only
fasta_files <- fasta_files[file.exists(fasta_files)]

if (length(fasta_files) > 0) {
  fasta <- readAAStringSet(fasta_files)
  
  # N-end rule amino acid categories
  type1     <- c("R", "K", "H")
  type2     <- c("F", "W", "Y", "L", "I")
  secondary <- c("D", "E")
  tertiary  <- c("N", "Q")
  small_res <- c("A", "C", "G", "P", "S", "T", "V")
  
  ntermini <- tibble(
    header   = names(fasta),
    sequence = as.character(fasta)
  ) %>%
    mutate(
      protein_length  = nchar(sequence),
      nterm_peptide   = substr(sequence, 1, 30),
      nterm_residue   = substr(sequence, 1, 1),
      second_residue  = substr(sequence, 2, 2),
      met_removed     = second_residue %in% small_res,
      predicted_nterm = ifelse(met_removed, second_residue, nterm_residue),
      
      protein_id      = str_extract(header, "^[^ ]+"),
      gene_symbol     = str_extract(header, "(?<=gene_symbol:)[^ ]+"),
      
      nend_class = case_when(
        predicted_nterm %in% type1 ~ "type1_destabilizing",
        predicted_nterm %in% type2 ~ "type2_destabilizing",
        predicted_nterm %in% secondary ~ "secondary_destabilizing",
        predicted_nterm %in% tertiary ~ "tertiary_destabilizing",
        TRUE ~ "stabilizing"
      )
    )
  
  write_csv(ntermini, here("data", "processed", "mouse_proteome_nterm_table.csv"))
  
  # Merge N-term annotations to Ub site sequences
  nterm_annotated <- ub_ss %>%
    left_join(ntermini, by = c("PEP.StrippedSequence" = "nterm_peptide"))
}

################################################################################
# 6. Linear Attenuation Modeling (Steatosis ~ Site vs Site + UBR1)
################################################################################
site_cols <- setdiff(ub_cols, c("Steatosis_perc_24", "Ubr1"))

results_att <- map_dfr(site_cols, function(site) {
  if (!site %in% colnames(ubsites)) return(NULL)
  
  df_sub <- ubsites[, c("Steatosis_perc_24", site, "Ubr1")]
  df_sub <- df_sub[complete.cases(df_sub), ]
  if (nrow(df_sub) < 10) return(NULL)
  
  # Wrap site names safely with backticks to support special characters
  f1 <- as.formula(paste0("Steatosis_perc_24 ~ `", site, "`"))
  f2 <- as.formula(paste0("Steatosis_perc_24 ~ `", site, "` + Ubr1"))
  
  m1 <- tryCatch(lm(f1, data = df_sub), error = function(e) NULL)
  m2 <- tryCatch(lm(f2, data = df_sub), error = function(e) NULL)
  
  if (is.null(m1) || is.null(m2)) return(NULL)
  
  beta1 <- coef(m1)[paste0("`", site, "`")]
  if (is.na(beta1)) beta1 <- coef(m1)[site]
  
  beta2 <- coef(m2)[paste0("`", site, "`")]
  if (is.na(beta2)) beta2 <- coef(m2)[site]
  
  p1 <- summary(m1)$coefficients[2, "Pr(>|t|)"]
  p2 <- summary(m2)$coefficients[2, "Pr(>|t|)"]
  
  attenuation <- ifelse(beta1 != 0, (beta1 - beta2) / beta1, NA)
  
  tibble(
    site = site,
    beta_site_only = beta1,
    beta_adjusted = beta2,
    p_site_only = p1,
    p_adjusted = p2,
    attenuation = attenuation
  )
}) %>%
  mutate(
    category = case_when(
      attenuation > 0.7 ~ "UBR1-mediated",
      attenuation > 0.3 ~ "partial mediation",
      TRUE ~ "independent"
    )
  )

################################################################################
# 7. Gene Set Enrichment Analysis (GSEA)
################################################################################
stats_df <- ubnorm_tau_sig %>%
  filter(!is.na(tau_ubr1) & !is.na(Gene)) %>%
  group_by(Gene) %>%
  summarize(tau_ubr1 = max(tau_ubr1), .groups = "drop") %>%
  arrange(desc(tau_ubr1))

stats <- setNames(stats_df$tau_ubr1, stats_df$Gene)

# Load MSigDB GO-BP gene sets
genesets_df <- msigdbr(species = "Mus musculus", category = "C5", subcollection = "BP")
genesets <- split(genesets_df$gene_symbol, genesets_df$gs_name)

# Execute GSEA
gsea_results <- fgsea(
  pathways  = genesets,
  stats     = stats,
  scoreType = "pos",
  minSize   = 10,
  maxSize   = 500,
  nproc     = 1
)

gsea_sig <- gsea_results %>%
  filter(pval <= 0.05) %>%
  arrange(desc(NES))

# Plot Top 10 Enriched Pathways
ggplot(head(gsea_sig, 10), aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top Enriched Gene Sets (UBR1 Tau Ranked)",
    x = "Gene Set (GO:BP)",
    y = "Normalized Enrichment Score (NES)"
  )