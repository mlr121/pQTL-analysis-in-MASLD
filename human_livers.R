library(readr)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsignif)

#load data ----
hl_proteome <- read.csv("rawdata/human/human_proteomics.csv", stringsAsFactors = F)

human_metadata <- read.csv("rawdata/human/human_metadata.csv", stringsAsFactors = F)
human_metadata$Run <- as.character(human_metadata$Run)
human_metadata$Steatosis.grade <- as.numeric(human_metadata$Steatosis.grade)

data_human <- hl_proteome %>%
  mutate(PG.ProteinGroups = NULL) %>% mutate(PG.FastaHeaders = NULL) 
data_human$PG.Genes <- make.unique(data_human$PG.Genes)
row.names(data_human) <- data_human$PG.Genes
data_human$PG.Genes <- NULL
data_human <- data.frame(t(data_human))

data_human <- data_human %>%
  rownames_to_column() %>%
  mutate(rowname = NULL) %>%
  rownames_to_column() %>%
  left_join(human_metadata, by = c("rowname" = "Run")) %>%
  mutate(meds = if_else(Diabetic.medications == "N/A", FALSE, TRUE))

#violin plot UBR1 expression ----
ggplot(data_human, aes(x = Diagnosis, y = log(UBR1))) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("NASH", " not NAFLD")), 
              map_signif_level=TRUE) +
  geom_point(aes(color = as.factor(NAS.score), shape = meds)) + 
  theme_minimal()

#tau ----
pheno_cols <- colnames(human_metadata)
pheno_cols <- pheno_cols[sapply(human_metadata[pheno_cols], is.numeric)]

compute_gam_nocov <- function(all_pheno, pheno_cols) {
  # Requires mgcv
  if (!"mgcv" %in% .packages(all.available = TRUE)) {
    suppressPackageStartupMessages(require(mgcv))
  } else {
    suppressPackageStartupMessages(library(mgcv))
  }
  
  # Basic checks
  if (!all(pheno_cols %in% colnames(all_pheno))) {
    missing_cols <- pheno_cols[!pheno_cols %in% colnames(all_pheno)]
    stop("These pheno_cols are not in all_pheno: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # Identify numeric columns
  numeric_cols <- colnames(all_pheno)[sapply(all_pheno, is.numeric)]
  
  # Proteins = numeric columns that are NOT phenotypes
  protein_cols <- setdiff(numeric_cols, pheno_cols)
  
  if (length(protein_cols) == 0)
    stop("No protein columns detected. Check pheno_cols selection.")
  
  results_list <- vector("list", length = length(pheno_cols))
  names(results_list) <- pheno_cols
  
  # Loop through phenotypes
  for (trait in pheno_cols) {
    message("Processing trait: ", trait)
    y <- all_pheno[[trait]]
    
    trait_res <- lapply(protein_cols, function(prot) {
      protein_vals <- all_pheno[[prot]]
      
      # Require >= 3 paired non-missing observations
      if (sum(!is.na(protein_vals) & !is.na(y)) < 3) return(NULL)
      
      df_model <- data.frame(
        trait = y,
        protein = protein_vals
      )
      
      fit <- tryCatch({
        mgcv::gam(trait ~ s(protein), data = df_model)
      }, error = function(e) return(NULL))
      
      if (is.null(fit)) return(NULL)
      
      s_summary <- summary(fit)
      
      pval <- tryCatch({
        s_summary$s.table["s(protein)", "p-value"]
      }, error = function(e) NA_real_)
      
      dev_expl <- tryCatch({
        s_summary$dev.expl
      }, error = function(e) NA_real_)
      
      tau_val <- tryCatch({
        cor(df_model$protein, df_model$trait,
            method = "kendall", use = "complete.obs")
      }, error = function(e) NA_real_)
      
      direction <- if (is.na(tau_val)) {
        NA_character_
      } else if (tau_val > 0) {
        "1"
      } else if (tau_val < 0) {
        "-1"
      } else {
        "neutral"
      }
      
      data.frame(
        Protein = prot,
        p_value = pval,
        deviance_explained = dev_expl,
        tau = tau_val,
        direction = direction,
        stringsAsFactors = FALSE
      )
    })
    
    trait_df <- do.call(rbind, trait_res)
    results_list[[trait]] <- if (!is.null(trait_df) && nrow(trait_df) > 0) trait_df else NULL
  }
  
  return(results_list)
}

##run ----
human_tau_allpheno <- compute_gam_nocov(data_human, pheno_cols)

#Combine
human_tau <- do.call(rbind, lapply(names(human_tau_allpheno), function(trait) {
  df <- human_tau_allpheno[[trait]]
  if (!is.null(df)) df$Trait <- trait
  df
}))

#mouse orthologs ----
library(biomaRt)

human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
human_genes <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol"
  ),
  mart = human_mart
)

colnames(human_genes) <- c(
  "human_ensembl",
  "human_symbol"
)
human_mouse_orth <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "mmusculus_homolog_ensembl_gene",
    "mmusculus_homolog_associated_gene_name",
    "mmusculus_homolog_orthology_type"
  ),
  mart = human_mart
)

colnames(human_mouse_orth) <- c(
  "human_ensembl",
  "mouse_ensembl",
  "mouse_symbol",
  "orthology_type"
)
orth_map <- human_genes |>
  dplyr::inner_join(human_mouse_orth, by = "human_ensembl")


#use if only 1-1 desired
#orth_1to1 <- orth_map %>%
#  dplyr::filter(orthology_type == "ortholog_one2one")

human_genes_present <- hl_proteome$PG.Genes
mouse_genes_present <- gsub("__\\d+", "", proteomics$Gene.Symbol)

shared_map <- orth_map %>%
  dplyr::filter(
    human_symbol %in% human_genes_present,
    mouse_symbol %in% mouse_genes_present)
shared_map[shared_map == ""] <- NA
shared_map <- shared_map %>%
  filter(!is.na(human_symbol))

human_tau_clean <- human_tau %>%
  dplyr::select(
    human_symbol = Protein,
    human_tau = tau,
    human_p = p_value,
    human_dev = deviance_explained,
    human_dir = direction,
    human_trait = Trait
  )

mouse_tau_clean <- combined_df %>%
  dplyr::select(
    mouse_symbol = Protein,
    mouse_tau = tau,
    mouse_p = p_value,
    mouse_dev = deviance_explained,
    mouse_dir = direction,
    mouse_trait = Trait
  )

##join mouse and human tau ----
master_tau <- orth_map %>%
  dplyr::inner_join(human_tau_clean, by = "human_symbol") %>%
  dplyr::inner_join(mouse_tau_clean, by = "mouse_symbol")

master_tau <- master_tau %>%
  dplyr::mutate(
    sign_concordant = dplyr::case_when(
      is.na(human_tau) | is.na(mouse_tau) ~ NA,
      human_tau * mouse_tau > 0 ~ TRUE,
      human_tau * mouse_tau < 0 ~ FALSE,
      TRUE ~ NA
    ),
    abs_tau_human = abs(human_tau),
    abs_tau_mouse = abs(mouse_tau)
  )

pqtl_clean <- pqtl_topLOD %>%
  dplyr::select(
    mouse_symbol = Phenotype,
    Chr,
    qtl_type,
    pos_bp = Pos_bp,
    Lod
  ) %>%
  dplyr::distinct(mouse_symbol, .keep_all = TRUE)

master_tau <- master_tau %>%
  dplyr::left_join(pqtl_clean, by = "mouse_symbol") %>%
  dplyr::mutate(
    has_pqtl = !is.na(Lod),
    cis_pqtl = qtl_type == "cis",
    trans_pqtl = qtl_type == "trans"
  )

##parse to one trait ----
steatosis_tau <- master_tau %>%
  dplyr::filter(
    human_trait == "Liver.fat....",
    mouse_trait == "Steatosis_perc_24"
  ) %>% filter(has_pqtl == TRUE)

##plot concordance ----

scatter <- ggplot(steatosis_tau,
                  aes(x = human_tau, y = mouse_tau, color = cis_pqtl)) +
  geom_point(aes(text = paste("Protein: ", human_symbol), size = Lod)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("lightgrey", "darkgreen")) +
  theme_minimal()
ggplotly(scatter)

#PTM search ----
human_ptm <- read.csv("rawdata/human/20260108_102414_HumanLiver_Report_PTM.csv", stringsAsFactors = F)

human_ub <- human_ptm %>%
  filter(grepl("GlyGly", EG.PrecursorId)) %>%
  filter(PEP.IsGeneSpecific == "TRUE") %>%
  mutate(EG.ProteinPTMLocations = NULL) %>%
  mutate(PEP.IsProteinGroupSpecific = NULL) %>%
  mutate(PEP.IsGeneSpecific = NULL) %>%
  mutate(PG.Genes = paste(PG.Genes, EG.PrecursorId, sep = "")) %>%
  mutate(EG.PrecursorId = NULL) %>%
  mutate(across(
    where(is.character),
    ~ gsub("Filtered", "", .x)
  ))

row.names(human_ub) <- human_ub$PG.Genes
human_ub$PG.Genes <- NULL
human_ub <- t(human_ub)
row.names(human_ub) <- gsub("X", "", row.names(human_ub))
human_ub <- data.frame(human_ub) %>%
  mutate(across(
    everything(),
    ~ as.numeric(.x)
  ))
#correlation to ubr1 abundance ----
human_ub_df <- data.frame(human_ub) %>%
  rownames_to_column() %>%
  left_join(data_human %>% dplyr::select(UBR1, rowname), by = "rowname")

pheno_cols <- "UBR1"
ubr1_ubsite_tau_human <- compute_gam_nocov(human_ub_df, pheno_cols = "UBR1")
ubr1_ubsite_tau_human <- as.data.frame(ubr1_ubsite_tau_human) %>%
  rename_with(~str_remove_all(., "UBR1."))
#ptm tau gsea ----

library(msigdbr)
library(fgsea)
library(tidyr)
ubr1_ubsite_tau_human <- ubr1_ubsite_tau_human %>% mutate(Gene = sub("_.*$", "", Protein))

data_ranked <- ubr1_ubsite_tau_human[order(ubr1_ubsite_tau_human$tau), ]

stats <- data_ranked$tau
names(stats) <- data_ranked$Gene  

# remove NA + duplicates
stats <- stats[!is.na(stats)]
stats <- tapply(stats, names(stats), max)

# sort decreasing
stats <- sort(stats, decreasing = TRUE)

# gene sets
geneset <- c(msigdbr(species = "human", category = "C5"
                   , subcollection = c("BP"
)), msigdbr(species = "human", category = "C5"
           , subcollection = c("MF"
           )))
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
print(gsea_sig$leadingEdge[order(gsea_sig$NES, decreasing = TRUE)])

#pathway enrichment
library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(dplyr)

# Gene sets (C5 BP + MF)
geneset <- rbind(
  msigdbr(species = "human", category = "C5", subcollection = "BP")
 # , msigdbr(species = "human", category = "C5", subcollection = "MF")
)
geneset_df <- geneset %>% dplyr::select(gs_name, gene_symbol)

# Run ORA enrichment
ora <- enricher(
  gene       = mouse_genes_present,   # your gene list
  TERM2GENE  = geneset_df,
  minGSSize  = 10,
  maxGSSize  = 500,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

ora_df <- as.data.frame(ora)

# Clean up pathway names for plotting
ora_df <- ora_df %>%
  mutate(
    pathway_clean = gsub("^GOBP_|^GOMF_", "", Description),
    pathway_clean = gsub("_", " ", pathway_clean),
    pathway_clean = stringr::str_to_sentence(pathway_clean)
  )

# Plot top 20
ora_df %>%
  slice_min(p.adjust, n = 20) %>%
  mutate(GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text = x)))) %>%
  ggplot(aes(x = reorder(pathway_clean, GeneRatio_num), 
             y = GeneRatio_num, 
             fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "#2166ac", high = "#d1e5f0",
                      name = "Adjusted\np-value") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.y   = element_text(size = 10),
    axis.text.x   = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Pathway Enrichment Analysis (ORA)",
    x     = NULL,
    y     = "Gene Ratio"
  )

#compare to protein level tau ----
ubr1_ub_tau_human <- ubr1_ubsite_tau_human %>%
  mutate(
    Gene = sub("_.*", "", Protein)
  ) %>%
  left_join(human_steatosis_tau %>% dplyr::select(Protein, tau), by = c("Gene" = "Protein")) %>%
  rename(ubr1_ub_tau = tau.x) %>%
  rename(stea_tau = tau.y)


ubr1_ub <- ggplot(ubr1_ub_tau_human, aes(x = ubr1_ub_tau, y = stea_tau)) +
  geom_point(aes(text = paste(Protein, Gene))) +
  theme_minimal()
ggplotly(ubr1_ub)


human_ub_df <- human_ub_df %>%
  left_join(data_human, by = "rowname")

ggplot(human_ub_df, aes(x = UBR1, y = PID1_PVIELWK.LeuArgGlyGly._.3)) +
  geom_point(aes(shape = meds)) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),  # cubic spline
    se = FALSE,                      # show confidence band
    aes(group = Sex)
  ) +
  theme_minimal()