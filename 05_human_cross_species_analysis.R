################################################################################
# Script: 05_human_cross_species_analysis.R
# Description: Human liver proteomics processing, UBR1 association analysis,
#              human-mouse ortholog cross-species concordance, and ORA/GSEA.
################################################################################

################################################################################
# 1. Reproducibility & Environment Setup
################################################################################
set.seed(123)
options(stringsAsFactors = FALSE)

# Core Libraries
library(readr)
library(tibble)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggsignif)
library(plotly)
library(here)

# Functional & Genomics Libraries
library(biomaRt)
library(msigdbr)
library(fgsea)
library(clusterProfiler)

# Source existing helper functions
source(here("functions", "gam_helpers.R"))

################################################################################
# 2. Ingest Human Proteomics & Metadata
################################################################################
hl_proteome <- read.csv(here("data", "raw", "human", "human_proteomics.csv"))
human_metadata <- read.csv(here("data", "raw", "human", "human_metadata.csv"))

human_metadata <- human_metadata %>%
  mutate(
    Run = as.character(Run),
    Steatosis.grade = as.numeric(Steatosis.grade)
  )

# Process proteomics matrix
data_human_mat <- hl_proteome %>%
  select(-any_of(c("PG.ProteinGroups", "PG.FastaHeaders"))) %>%
  mutate(PG.Genes = make.unique(PG.Genes)) %>%
  column_to_rownames("PG.Genes") %>%
  t() %>%
  as.data.frame()

data_human <- data_human_mat %>%
  rownames_to_column(var = "Run") %>%
  left_join(human_metadata, by = "Run") %>%
  mutate(
    meds = if_else(Diabetic.medications == "N/A" | is.na(Diabetic.medications), FALSE, TRUE)
  )

################################################################################
# 3. Exploratory UBR1 Expression Visualization
################################################################################
p_ubr1_violin <- ggplot(data_human, aes(x = Diagnosis, y = log(UBR1))) +
  geom_boxplot(outlier.shape = NA) +
  geom_signif(
    comparisons = list(c("NASH", " not NAFLD")),
    map_signif_level = TRUE
  ) +
  geom_point(aes(color = as.factor(NAS.score), shape = meds), position = position_jitter(width = 0.15)) +
  theme_minimal() +
  labs(title = "Human UBR1 Expression by Diagnosis", color = "NAS Score", shape = "On Meds")

print(p_ubr1_violin)

################################################################################
# 4. Human Phenotype Associations (GAM & Kendall's Tau)
################################################################################
pheno_cols <- colnames(human_metadata)[sapply(human_metadata, is.numeric)]

human_tau_allpheno <- compute_gam_nocov(data_human, pheno_cols = pheno_cols)

# Flatten list to single data frame
human_tau <- bind_rows(lapply(names(human_tau_allpheno), function(trait) {
  df <- human_tau_allpheno[[trait]]
  if (!is.null(df)) df$Trait <- trait
  df
}))

################################################################################
# 5. Human-Mouse Ortholog Mapping via biomaRt
################################################################################
tryCatch({
  human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  human_genes <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    mart = human_mart
  ) %>%
    rename(human_ensembl = ensembl_gene_id, human_symbol = hgnc_symbol)

  human_mouse_orth <- getBM(
    attributes = c(
      "ensembl_gene_id",
      "mmusculus_homolog_ensembl_gene",
      "mmusculus_homolog_associated_gene_name",
      "mmusculus_homolog_orthology_type"
    ),
    mart = human_mart
  ) %>%
    rename(
      human_ensembl = ensembl_gene_id,
      mouse_ensembl = mmusculus_homolog_ensembl_gene,
      mouse_symbol = mmusculus_homolog_associated_gene_name,
      orthology_type = mmusculus_homolog_orthology_type
    )

  orth_map <- human_genes %>%
    inner_join(human_mouse_orth, by = "human_ensembl") %>%
    filter(!is.na(human_symbol) & human_symbol != "")
    
  write.csv(orth_map, here("data", "processed", "human_mouse_ortholog_map.csv"), row.names = FALSE)
}, error = function(e) {
  message("⚠️ biomaRt query failed or offline. Attempting to load cached map.")
  if (file.exists(here("data", "processed", "human_mouse_ortholog_map.csv"))) {
    orth_map <<- read.csv(here("data", "processed", "human_mouse_ortholog_map.csv"))
  }
})

################################################################################
# 6. Cross-Species Concordance Analysis
################################################################################
if (exists("orth_map") && exists("combined_df") && exists("pqtl_topLOD")) {
  
  human_tau_clean <- human_tau %>%
    select(
      human_symbol = Protein,
      human_tau = tau,
      human_p = p_value,
      human_dev = deviance_explained,
      human_dir = direction,
      human_trait = Trait
    )

  mouse_tau_clean <- combined_df %>%
    select(
      mouse_symbol = Protein,
      mouse_tau = tau,
      mouse_p = p_value,
      mouse_dev = deviance_explained,
      mouse_dir = direction,
      mouse_trait = Trait
    )

  master_tau <- orth_map %>%
    inner_join(human_tau_clean, by = "human_symbol") %>%
    inner_join(mouse_tau_clean, by = "mouse_symbol") %>%
    mutate(
      sign_concordant = case_when(
        is.na(human_tau) | is.na(mouse_tau) ~ NA,
        human_tau * mouse_tau > 0 ~ TRUE,
        human_tau * mouse_tau < 0 ~ FALSE,
        TRUE ~ NA
      ),
      abs_tau_human = abs(human_tau),
      abs_tau_mouse = abs(mouse_tau)
    )

  pqtl_clean <- pqtl_topLOD %>%
    select(mouse_symbol = Phenotype, Chr, qtl_type, pos_bp = Pos_bp, Lod) %>%
    distinct(mouse_symbol, .keep_all = TRUE)

  master_tau <- master_tau %>%
    left_join(pqtl_clean, by = "mouse_symbol") %>%
    mutate(
      has_pqtl = !is.na(Lod),
      cis_pqtl = qtl_type == "cis",
      trans_pqtl = qtl_type == "trans"
    )

  steatosis_tau <- master_tau %>%
    filter(
      human_trait == "Liver.fat....",
      mouse_trait == "Steatosis_perc_24",
      has_pqtl == TRUE
    )

  # Concordance Scatter Plot
  p_scatter <- ggplot(steatosis_tau, aes(x = human_tau, y = mouse_tau, color = cis_pqtl)) +
    geom_point(aes(text = paste("Protein:", human_symbol), size = Lod)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "lightgrey")) +
    theme_minimal() +
    labs(x = "Human Tau (Liver Fat)", y = "Mouse Tau (Steatosis 24wk)")

  ggplotly(p_scatter)
}

################################################################################
# 7. Human PTM (GlyGly / Ubiquitination) Analysis
################################################################################
human_ptm <- read.csv(here("data", "raw", "human", "20260108_102414_HumanLiver_Report_PTM.csv"))

human_ub <- human_ptm %>%
  filter(str_detect(EG.PrecursorId, "GlyGly")) %>%
  filter(PEP.IsGeneSpecific == TRUE | PEP.IsGeneSpecific == "TRUE") %>%
  mutate(
    PG.Genes = paste0(PG.Genes, EG.PrecursorId),
    across(where(is.character), ~ str_remove_all(., "Filtered"))
  ) %>%
  select(-c(EG.ProteinPTMLocations, PEP.IsProteinGroupSpecific, PEP.IsGeneSpecific, EG.PrecursorId)) %>%
  column_to_rownames("PG.Genes") %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric))

rownames(human_ub) <- str_remove_all(rownames(human_ub), "^X")

human_ub_df <- human_ub %>%
  rownames_to_column(var = "Run") %>%
  left_join(data_human %>% select(UBR1, Run, meds, Sex), by = "Run")

# UBR1 GAM associations against PTM sites
ubr1_ubsite_tau_human <- compute_gam_nocov(human_ub_df, pheno_cols = "UBR1")
ubr1_ubsite_tau_human <- as.data.frame(ubr1_ubsite_tau_human$UBR1)

################################################################################
# 8. Human PTM Functional Enrichment (GSEA & ORA)
################################################################################
# GSEA Pipeline
ubr1_ubsite_tau_human <- ubr1_ubsite_tau_human %>%
  mutate(Gene = sub("_.*$", "", Protein))

stats_df <- ubr1_ubsite_tau_human %>%
  filter(!is.na(tau) & !is.na(Gene)) %>%
  group_by(Gene) %>%
  summarize(tau = max(tau), .groups = "drop") %>%
  arrange(desc(tau))

stats <- setNames(stats_df$tau, stats_df$Gene)

# Load combined BP + MF gene sets cleanly
geneset_df <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  filter(subcollection %in% c("BP", "MF"))

geneset_list <- split(geneset_df$gene_symbol, geneset_df$gs_name)

gsea_res <- fgsea(
  pathways  = geneset_list,
  stats     = stats,
  scoreType = "pos",
  minSize   = 10,
  maxSize   = 500,
  nproc     = 1
)

gsea_sig <- gsea_res %>%
  filter(pval <= 0.05) %>%
  arrange(desc(NES))

# Plot Top GSEA Pathways
ggplot(head(gsea_sig, 10), aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(fill = "darkred") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Human Ub-Site GSEA (UBR1 Tau)", x = "Pathway", y = "NES")

# Over-Representation Analysis (ORA) with clusterProfiler
if (exists("mouse_genes_present")) {
  ora_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcollection = "BP") %>%
    select(gs_name, gene_symbol)

  ora <- enricher(
    gene         = mouse_genes_present,
    TERM2GENE    = ora_sets,
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  if (!is.null(ora)) {
    ora_df <- as.data.frame(ora) %>%
      mutate(
        pathway_clean = str_to_sentence(gsub("_", " ", gsub("^GOBP_|^GOMF_", "", Description))),
        # Convert "x/y" ratio safely without eval(parse())
        GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
      )

    p_ora <- ora_df %>%
      slice_min(p.adjust, n = 20) %>%
      ggplot(aes(x = reorder(pathway_clean, GeneRatio_num), y = GeneRatio_num, fill = p.adjust)) +
      geom_col(width = 0.7) +
      scale_fill_gradient(low = "#2166ac", high = "#d1e5f0", name = "Adjusted p-value") +
      coord_flip() +
      theme_minimal(base_size = 12) +
      labs(title = "Pathway Enrichment Analysis (ORA)", x = NULL, y = "Gene Ratio")

    print(p_ora)
  }
}