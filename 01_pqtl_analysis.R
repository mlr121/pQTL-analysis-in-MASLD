################################################################################
# Script: 01_pqtl_analysis.R
# Description: Female pQTL processing, cis/trans annotation, plotting, and GAMs.
################################################################################

################################################################################
# 1. Reproducibility & Environment Setup
################################################################################
set.seed(123)
options(stringsAsFactors = FALSE)

# Load core libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(plotly)
library(patchwork)
library(here)
library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(BioVenn)
library(ggpmisc)

# Source helper functions
source(here("functions", "gam_helpers.R"))

################################################################################
# 2. pQTL Loading & Filtering
################################################################################
fdata <- read.delim(here("data", "raw", "final_Gene_Female_log2_filter05_ProMatrix_QTL.tsv"))

pqtls <- fdata %>%
  select(-Lod_threshold) %>%
  filter(Status %in% c("Significant_001", "Significant_005")) %>%
  mutate(Phenotype = str_remove_all(Phenotype, "_[2-5]"))

################################################################################
# 3. Cis/Trans Annotation
################################################################################
ensdb <- EnsDb.Mmusculus.v79
protein_locations <- as.data.frame(genes(ensdb)) %>%
  select(symbol, Gene_chr, start, end, strand)

# Standardize Chromosome Order
Chromosome_order <- c(as.character(1:19), "X", "Y")
protein_locations <- protein_locations %>%
  mutate(Gene_chr = as.integer(factor(Gene_chr, levels = Chromosome_order)))

# Merge and Annotate
cis_distance_threshold <- 1000000

pqtls_annotated <- pqtls %>%
  left_join(protein_locations, by = c("Phenotype" = "symbol")) %>%
  rename(Gene_start = start) %>% # renaming for clarity based on original script usage
  mutate(
    qtl_type = case_when(
      (Gene_chr == Chr) & (abs(Gene_start - Pos_bp) <= cis_distance_threshold) ~ "cis",
      (Gene_chr != Chr) | (abs(Gene_start - Pos_bp) > cis_distance_threshold) ~ "trans",
      TRUE ~ NA_character_
    )
  )

################################################################################
# 4. Top LOD Filtering & NASH Annotation
################################################################################
# Replaced data.table logic with dplyr::slice_max for consistency
pqtl_topLOD <- pqtls_annotated %>%
  group_by(Phenotype, Chr, qtl_type) %>%
  slice_max(order_by = Lod, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(Lod))

qtl_counts <- pqtl_topLOD %>% count(qtl_type)
print(qtl_counts)

nash_list <- tolower(c(
  "Pnpla3", "Tm6sf2", "Gckr", "Mboat7", "Hsd17b13__1", "Fasn", "Acaca", "Acly",
  "Scd1", "Srebf1", "Mlxipl", "Hmgcs1", "Hmgcr", "Dgat2", "Nr1h4", "Cd36",
  "Fabp1", "Fabp4", "Lpl", "Ldlr", "Plin2", "Apoe", "Apob", "Cpt1a", "Acox1",
  "Cyp2e1", "Ppara", "Pparg", "Thrb", "Hnf4a", "Pck1", "Col1a1", "Acta2",
  "Tgfb1", "Timp1", "Mmp2", "Mmp9", "Tnf", "Trem2", "Socs3", "Il10",
  "Fgf21", "Gdf15", "Adipoq", "Saa1", "Sod2", "Gpx1", "Hmox1", "Gpam", "Pnpla2"
))

pqtl_topLOD <- pqtl_topLOD %>%
  mutate(nash_protein = as.factor(ifelse(tolower(Phenotype) %in% nash_list, Phenotype, NA)))

################################################################################
# 5. Manhattan & Density Plotting
################################################################################
chr_lengths <- c(
  "1" = 195471971, "2" = 182113224, "3" = 160039680, "4" = 156508116, 
  "5" = 151834684, "6" = 149736546, "7" = 145441459, "8" = 129401213, 
  "9" = 124595110, "10"= 130694993, "11"= 122082543, "12"= 120129022,
  "13"= 120421639, "14"= 124902244, "15"= 104043685, "16"= 98207768,  
  "17"= 94987271,  "18"= 90702639,  "19"= 61431566
)
gap <- 0 
chr_offsets <- cumsum(c(0, chr_lengths[-length(chr_lengths)] + gap))
names(chr_offsets) <- names(chr_lengths)
chr_mids <- chr_offsets + chr_lengths / 2

chr_bounds <- tibble(
  Chr = names(chr_lengths),
  xmin = chr_offsets,
  xmax = chr_offsets + chr_lengths,
  shade = as.integer(names(chr_lengths)) %% 2 == 0
)

plot_df <- pqtl_topLOD %>%
  mutate(
    Chr = as.character(Chr),
    Gene_chr = as.character(Gene_chr),
    abs_pos = Pos_bp + chr_offsets[Chr],
    abs_gene = Gene_start + chr_offsets[Gene_chr]
  )

x_scale <- scale_x_continuous(breaks = chr_mids, labels = names(chr_mids), expand = c(0.01, 0))

p_density <- ggplot() +
  geom_density(data = filter(plot_df, qtl_type == "trans"), aes(x = abs_pos), color = "firebrick", linewidth = 0.8, adjust = 0.25) +
  geom_density(data = filter(plot_df, qtl_type == "cis"), aes(x = abs_pos), color = "steelblue", linewidth = 0.8, adjust = 0.25) +
  x_scale +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(y = "Density")

p_map <- ggplot() +
  geom_rect(data = chr_bounds, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = shade), alpha = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "grey90", "FALSE" = "white")) +
  geom_point(data = plot_df, aes(x = abs_pos, y = abs_gene, color = qtl_type), size = 2, alpha = 0.7) +
  x_scale +
  scale_y_continuous(breaks = chr_mids, labels = names(chr_mids), expand = c(0.01, 0)) +
  scale_color_manual(values = c("cis" = "steelblue", "trans" = "firebrick"), name = "QTL Type") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) +
  labs(x = "QTL Position", y = "Gene Position")

p_density / p_map + plot_layout(heights = c(1, 4))
# ggsave(here("output", "figures", "density_map.png"), width = 8, height = 6)

################################################################################
# 6. Male vs Female Comparison
################################################################################
# WARNING: mpqtl_topLOD is missing from the global environment in this snippet.
# You will need to load it here before running the next lines.
# mpqtl_topLOD <- read.csv(here("data", "processed", "male_pqtl_topLOD.csv"))

# male_selected <- mpqtl_topLOD[, c("Phenotype", "Chr", "qtl_type", "Lod")]
# female_selected <- pqtl_topLOD[, c("Phenotype", "Chr", "qtl_type", "Lod", "nash_protein")]
# mvf <- merge(male_selected, female_selected, by = c("Phenotype", "Chr", "qtl_type"), suffixes = c("_male", "_female"))

# male_data <- apply(mpqtl_topLOD[, c("Phenotype", "Chr", "qtl_type")], 1, paste, collapse = "_")
# female_data <- apply(pqtl_topLOD[, c("Phenotype", "Chr", "qtl_type")], 1, paste, collapse = "_")
# draw.venn(male_data, female_data, NULL, title = "shared QTLs", subtitle = NA, xtitle = "male", ytitle = "female")

################################################################################
# 7. Proteomics Data Prep & GAM Modeling
################################################################################
proteomics <- read.csv(here("data", "raw", "fixcol_ProMatrix_metadata.csv"))
row.names(proteomics) <- proteomics$Gene.Symbol
proteomics$Gene.Symbol <- NULL

protein_matrix <- as.data.frame(t(proteomics)) %>% mutate(ID = row.names(.))

pheno <- read.csv(here("data", "raw", "F2_phenotypes.csv")) %>%
  mutate(ID = str_replace_all(ID, "-", "_"))

non_traits <- c("ID", "Strain", "Generation", "Sex")
pheno_cols <- colnames(pheno)[!(colnames(pheno) %in% non_traits)]
pheno_cols <- pheno_cols[sapply(pheno[pheno_cols], is.numeric)]

all_pheno <- protein_matrix %>% left_join(pheno, by = "ID")

# Run GAM calculations (using the helper functions sourced at the top)
results_tau_allpheno <- compute_gam_associations(all_pheno, pheno_cols)

combined_df <- do.call(rbind, lapply(names(results_tau_allpheno), function(trait) {
  df <- results_tau_allpheno[[trait]]
  if (!is.null(df)) df$Trait <- trait
  df
}))

ubr1_results <- compute_ubr1_associations(all_pheno, ubr1_col = "Ubr1")