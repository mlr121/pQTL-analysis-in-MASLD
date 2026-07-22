# PTM search ----
##load data ----
ptm <- read.csv("rawdata/simple_UbReport.txt", stringsAsFactors = F, sep = "\t")

ub <- ptm %>%
  mutate(PEP.IsProteinGroupSpecific = NULL) %>%
  mutate(PG.Genes = paste(Gene, EG.PrecursorId, sep = "")) %>%
  mutate(EG.PrecursorId = NULL) %>%
  mutate(PG.FastaHeaders = NULL) %>%
  mutate(EG.ProteinPTMLocations = NULL) %>%
  mutate(Gene = NULL) %>%
  mutate(PEP.StrippedSequence = NULL)

row.names(ub) <- ub$PG.Genes
ub$PG.Genes <- NULL
ub <- t(ub)
row.names(ub) <- gsub("X", "", row.names(ub))
ub <- data.frame(ub) %>%
  mutate(across(
    everything(),
    ~ as.numeric(.x)
  ))

ub_ss <- ptm %>% select(c(Gene, EG.PrecursorId, PEP.StrippedSequence)) %>%
  mutate(PG.Genes = paste(Gene, EG.PrecursorId, sep = ""))


# load proteomics----
proteomics_raw <- read.csv(file="raw_data/fixcol_ProMatrix_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
pheno <- read.csv("rawdata/F2_phenotypes.csv", stringsAsFactors = F)

proteomics_raw$Gene.Symbol <- ifelse(
  is.na(proteomics_raw$Gene.Symbol) | proteomics_raw$Gene.Symbol == "",
  proteomics_raw$protein.groups,
  proteomics_raw$Gene.Symbol
)
# Remove columns not containing numeric abundance data
proteomics_clean <- dplyr::select(proteomics_raw, -c(Gene.Symbol, protein.groups, PG.FastaHeaders))

# Transpose: rows = samples, columns = proteins
proteomics <- t(proteomics_clean)
proteomics <- as.data.frame(proteomics)
rownames(proteomics) <- gsub("_", "", rownames(proteomics))

# Assign protein names from the original data
colnames(proteomics) <- make.unique(proteomics_raw$Gene.Symbol)

# Add sample IDs
proteomics$ID <- rownames(proteomics)
colnames(proteomics) <- make.unique(colnames(proteomics))

pheno <- pheno %>%
  mutate(ID = gsub("-", "", ID))

non_traits <- c("ID","Strain","Generation","Sex")

pheno_cols <- colnames(pheno)[!(colnames(pheno) %in% non_traits)]
pheno_cols <- pheno_cols[sapply(pheno[pheno_cols], is.numeric)]

proteomics$ID <- make.unique(row.names(proteomics))
colnames(proteomics) <- make.unique(colnames(proteomics))
all_pheno <- pheno %>%
  left_join(proteomics, by = "ID")


# stitch to all_pheno----
ub$ID <- row.names(ub)
ubsites <- ub %>%
  left_join(all_pheno, by = "ID")


# normalize ub sites by protein abundance ----
ub_log <- ubsites %>%
  dplyr::mutate(
    dplyr::across(
      .cols = where(is.numeric) & !all_of(pheno_cols),
      .fns  = log2
    )
  )

cn <- colnames(ub_log)
expr_cols <- setdiff(cn, pheno_cols)
ub_cols <- expr_cols[grepl("_", expr_cols)]
protein_cols <- setdiff(expr_cols, ub_cols)

ub_genes <- sub("_.*$", "", ub_cols)

has_protein <- ub_genes %in% protein_cols

ub_cols_matched  <- ub_cols[has_protein]
ub_genes_matched <- ub_genes[has_protein]

ub_norm <- ub_log

for (i in seq_along(ub_cols_matched)) {
  site_col <- ub_cols_matched[i]
  prot_col <- ub_genes_matched[i]
  
  ub_norm[[site_col]] <-
    ub_log[[site_col]] - ub_log[[prot_col]]
}



# plotting scatter ----
protein <- ggplot(all_pheno, aes(x = Ubr1, y = Plin2, alpha = Steatosis_perc_24)) +
  geom_point(aes(color = Sex)) +
  geom_smooth() +
  theme_minimal()

site <- ggplot(ub_norm, aes(x = Ubr1, y = 	
                              Mcmbp_LQHINPLLPTC.Carbamidomethyl..C..LNK.GlyGly..K..EESR_.4, alpha = Steatosis_perc_24)) +
  geom_point(aes(color = Sex), size = 3) +
  theme_minimal()
stea <- ggplot(ubsites, aes(x = Steatosis_perc_24, y = PID1)) +
  geom_point(aes(color = Sex)) +
  theme_minimal()

ggplotly(protein) 
ggplotly(site)
stea

# ubr1 correlations ----
compute_ubr1_associations <- function(all_pheno, ubr1_col = "Ubr1") {
  
  # Requires mgcv
  if (!"mgcv" %in% .packages(all.available = TRUE)) {
    suppressPackageStartupMessages(require(mgcv))
   } else {
    suppressPackageStartupMessages(library(mgcv))
   }
  
  # Basic checks
  if (!ubr1_col %in% colnames(all_pheno)) {
    stop("ubr1_col not found in all_pheno: ", ubr1_col)
   }
  if (!"Sex" %in% colnames(all_pheno)) {
    stop("all_pheno must contain a 'Sex' column.")
   }
  
  sex_factor <- as.factor(all_pheno$Sex)
  
  # Identify numeric columns
  numeric_cols <- colnames(all_pheno)[sapply(all_pheno, is.numeric)]
  
  # Outcomes = all numeric proteins except Ubr1
  protein_cols <- setdiff(numeric_cols, ubr1_col)
  
  if (length(protein_cols) == 0) {
    stop("No protein outcome columns detected.")
   }
  
  ubr1_vals <- all_pheno[[ubr1_col]]
  
  message("Computing associations for ", ubr1_col)
  
  results <- lapply(protein_cols, function(prot) {
    
    y <- all_pheno[[prot]]
    
    # Require >= 3 paired non-missing observations
    if (sum(!is.na(ubr1_vals) & !is.na(y)) < 3) return(NULL)
    
    df_model <- data.frame(
      outcome = y,
      sex = sex_factor,
      ubr1 = ubr1_vals
    )
    
    fit <- tryCatch({
      mgcv::gam(outcome ~ s(ubr1) + sex, data = df_model)
     }, error = function(e) return(NULL))
    
    if (is.null(fit)) return(NULL)
    
    s_summary <- summary(fit)
    
    pval <- tryCatch({
      s_summary$s.table["s(ubr1)", "p-value"]
     }, error = function(e) NA_real_)
    
    dev_expl <- tryCatch({
      s_summary$dev.expl
     }, error = function(e) NA_real_)
    
    tau_val <- tryCatch({
      cor(df_model$ubr1, df_model$outcome,
          method = "kendall", use = "complete.obs")
     }, error = function(e) NA_real_)
    
    direction <- if (is.na(tau_val)) NA_character_
    else if (tau_val > 0) "1"
    else if (tau_val < 0) "-1"
    else "neutral"
    
    data.frame(
      Protein = prot,
      p_value = pval,
      deviance_explained = dev_expl,
      tau = tau_val,
      direction = direction,
      stringsAsFactors = FALSE
    )
   })
  
  # bind + drop NULLs
  results_df <- do.call(rbind, results)
  
  return(results_df)
    }

## run----
ubr1_results <- compute_ubr1_associations(ub_norm, ubr1_col = "Ubr1")
stea_tau_results <- compute_ubr1_associations(ub_norm, ubr1_col = "Steatosis_perc_24")

#before norm
ubr1_results_nonorm <- compute_ubr1_associations(ubsites, ubr1_col = "Ubr1")
stea_tau_results_nonorm <- compute_ubr1_associations(ubsites, ubr1_col = "Steatosis_perc_24")


# bind to protein level tau----
ubnorm_tau <- stea_tau_results %>%
  filter(Protein %in% ub_cols) %>%
  mutate(
    Gene = sub("_.*", "", Protein)
  ) %>% left_join(ubr1_results, by = "Protein") %>%
  rename(tau_stea = tau.x) %>% rename(tau_ubr1 = tau.y)
ubnorm_tau_sig <- ubnorm_tau %>%
  dplyr::filter(p_value.y <= 0.05)

ubtau <- ggplot(ubnorm_tau, aes(x = tau_stea, y = tau_ubr1)) +
  geom_point(aes(text = paste(Protein))) +
  theme_minimal()

ggplotly(ubtau)


# bind to pQTL----
pqtl_topLOD_h <- pqtl_topLOD %>%
  left_join(steatosis_tau %>% select(human_tau, mouse_symbol, mouse_tau, sign_concordant), by = c("Phenotype" = "mouse_symbol")) %>%
  filter(!is.na(human_tau)) %>%
  filter(sign_concordant == TRUE)

h_qtl_tau <- ggplot(pqtl_topLOD_h, aes(x = qtl_loc, y = Lod)) +
  geom_point(aes(alpha = abs(human_tau), color = qtl_type, size = human_tau, text = paste(Phenotype))) +
  theme_minimal()

ggplotly(h_qtl_tau)  

# search for n termini ----
library(Biostrings)
library(tidyverse)

# load FASTA
fasta_file <- c("Mus_musculus_c57bl6nj.C57BL_6NJ_v1.pep.all.fasta", "Mus_musculus_casteij.CAST_EiJ_v1.pep.all.fasta", "Mus_musculus_pwkphj.PWK_PhJ_v1.pep.all.fasta", "Mus_musculus_129s1svimj.129S1_SvImJ_v1.pep.all.fasta")

# Length of N-terminal peptide to extract
n_len <- 30

# Read FASTA
fasta <- readAAStringSet(fasta_file)

# N-end rule groups
type1 <- c("R","K","H")
type2 <- c("F","W","Y","L","I")
secondary <- c("D","E")
tertiary <- c("N","Q")
small_res <- c("A","C","G","P","S","T","V")

ntermini <- tibble(
  header = names(fasta),
  sequence = as.character(fasta)
) %>%
  mutate(
    protein_length = nchar(sequence),
    nterm_peptide = substr(sequence, 1, n_len),
    nterm_residue = substr(sequence, 1, 1),
    second_residue = substr(sequence, 2, 2),
    
    # Predict initiator Met removal
    met_removed = second_residue %in% small_res,
    predicted_nterm = ifelse(met_removed,
                             substr(sequence, 2, 2),
                             substr(sequence, 1, 1)),
    
    # Parse Ensembl identifiers
    protein_id = str_extract(header, "^[^ ]+"),
    gene_id = str_extract(header, "(?<=gene:)[^ ]+"),
    transcript_id = str_extract(header, "(?<=transcript:)[^ ]+"),
    gene_symbol = str_extract(header, "(?<=gene_symbol:)[^ ]+"),
    
    # N-end rule classification
    nend_class = case_when(
      predicted_nterm %in% type1 ~ "type1_destabilizing",
      predicted_nterm %in% type2 ~ "type2_destabilizing",
      predicted_nterm %in% secondary ~ "secondary_destabilizing",
      predicted_nterm %in% tertiary ~ "tertiary_destabilizing",
      TRUE ~ "stabilizing"
    )
  ) %>%
  select(
    protein_id,
    gene_id,
    transcript_id,
    gene_symbol,
    protein_length,
    nterm_peptide,
    predicted_nterm,
    nend_class
  )

# Save output
write_csv(ntermini, "mouse_proteome_nterm_table.csv")

# Quick summary
table(ntermini$nend_class)


# stripped seq filtered by nterm ---- 
nterm <- ub_ss %>%
  left_join(ntermini, by = c("PEP.StrippedSequence" = "nterm_peptide"))


# linear model ----
model <- lm(Steatosis_perc_24 ~ Hsd17b13_RGVEETADK.GlyGly..K..C.Carbamidomethyl..C..R_.3 + Ubr1, data = ubsites)
summary(model)
model2 <- lm(Steatosis_perc_24 ~ Hsd17b13_RGVEETADK.GlyGly..K..C.Carbamidomethyl..C..R_.3, data = ubsites)
summary(model2)

# attenuation ----
library(broom)
library(dplyr)
library(purrr)

site_cols <- setdiff(colnames(ubsites), c("Steatosis_perc_24", "Ubr1"))

results_att <- map_dfr(site_cols, function(site) {
  
  # build formulas
  f1 <- as.formula(paste("Steatosis_perc_24 ~", site))
  f2 <- as.formula(paste("Steatosis_perc_24 ~", site, "+ Ubr1"))
  
  # fit models
  m1 <- lm(f1, data = ubsites)
  m2 <- lm(f2, data = ubsites)
  
  # extract coefficients
  beta1 <- coef(m1)[site]
  beta2 <- coef(m2)[site]
  
  p1 <- summary(m1)$coefficients[site, "Pr(>|t|)"]
  p2 <- summary(m2)$coefficients[site, "Pr(>|t|)"]
  
  # attenuation
  attenuation <- (beta1 - beta2) / beta1
  
  tibble(
    site = site,
    beta_site_only = beta1,
    beta_adjusted = beta2,
    p_site_only = p1,
    p_adjusted = p2,
    attenuation = attenuation
  )
 })
results_att <- results_att %>%
  mutate(
    category = case_when(
      attenuation > 0.7 ~ "UBR1-mediated",
      attenuation > 0.3 ~ "partial mediation",
      TRUE ~ "independent"
    )
  )


# chr2 qtls ----
chr2qtl <- pqtl_topLOD %>%
  dplyr::filter(chr == 2)
mchr2qtl <- mpqtl_topLOD %>%
  dplyr::filter(chr == 2)

chr2_ubtau <- ubnorm_tau %>%
  dplyr::filter(Gene %in% c(chr2qtl$Phenotype, mchr2qtl$Phenotype))

count(chr2_ubtau %>% dplyr::filter(tau_ubr1 > 0.5))
count(ubnorm_tau %>% dplyr::filter(tau_ubr1 > 0.5))

print(chr2_ubtau %>% dplyr::filter(tau_ubr1 > 0.5))

# gsea ----

library(msigdbr)
library(fgsea)
library(tidyr)

data_ranked <- ubnorm_tau_sig[order(-ubnorm_tau_sig$tau_ubr1), ]

stats <- data_ranked$tau_ubr1
names(stats) <- data_ranked$Gene  

# remove NA + duplicates
stats <- stats[!is.na(stats)]
stats <- tapply(stats, names(stats), max)

# sort decreasing
stats <- sort(stats, decreasing = TRUE)

# gene sets
geneset <- msigdbr(species = "Mus musculus", category = "C5", subcollection = "BP")
geneset <- split(geneset$gene_symbol, geneset$gs_name)

# run GSEA
GSEA <- fgsea(
  pathways = geneset, 
  stats = stats,
  scoreType = "pos", 
  minSize = 10,
  maxSize = 500,
  nproc = 1
)

gsea_sig <- GSEA %>%
  dplyr::filter(pval <= 0.05)

ggplot(gsea_sig[order(gsea_sig$NES, decreasing = TRUE), ][1:10, ], 
       aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Enriched Gene Sets", x = "Gene Set", y = "Normalized Enrichment Score (NES)")
