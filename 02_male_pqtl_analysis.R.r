################################################################################
# Script: 02_male_pqtl_analysis.R
# Description: Male pQTL processing, cis/trans annotation, and plotting.
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

################################################################################
# 2. pQTL Loading & Filtering
################################################################################
mdata <- read.delim(here("data", "raw", "final_Gene_Male_log2_filter05_ProMatrix_QTL.tsv"))

pqtls <- mdata %>%
  select(-Lod_threshold) %>%
  filter(Status %in% c("Significant_001", "Significant_005")) %>%
  mutate(Phenotype = str_remove_all(Phenotype, "_[2-5]")) # Cleaning suffixes if present

################################################################################
# 3. Genomic Annotation & Cis/Trans Classification
################################################################################
# Extract protein locations
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
  rename(Gene_start = start) %>%
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
# Filter to max LOD per phenotype/chromosome/type
mpqtl_topLOD <- pqtls_annotated %>%
  group_by(Phenotype, Chr, qtl_type) %>%
  slice_max(order_by = Lod, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(Lod))

# Print summary
qtl_counts <- mpqtl_topLOD %>% count(qtl_type)
print(qtl_counts)

# Annotate NASH proteins
nash_list <- tolower(c(
  "Pnpla3", "Tm6sf2", "Gckr", "Mboat7", "Hsd17b13__1", "Fasn", "Acaca", "Acly",
  "Scd1", "Srebf1", "Mlxipl", "Hmgcs1", "Hmgcr", "Dgat2", "Nr1h4", "Cd36",
  "Fabp1", "Fabp4", "Lpl", "Ldlr", "Plin2", "Apoe", "Apob", "Cpt1a", "Acox1",
  "Cyp2e1", "Ppara", "Pparg", "Thrb", "Hnf4a", "Pck1", "Col1a1", "Acta2",
  "Tgfb1", "Timp1", "Mmp2", "Mmp9", "Tnf", "Trem2", "Socs3", "Il10",
  "Fgf21", "Gdf15", "Adipoq", "Saa1", "Sod2", "Gpx1", "Hmox1", "Gpam", "Pnpla2"
))

mpqtl_topLOD <- mpqtl_topLOD %>%
  mutate(nash_protein = as.factor(ifelse(tolower(Phenotype) %in% nash_list, Phenotype, NA)))

# Optional: Save intermediate data for cross-comparison in script 01
# write.csv(mpqtl_topLOD, here("data", "processed", "male_pqtl_topLOD.csv"), row.names = FALSE)

################################################################################
# 5. Manhattan Plot
################################################################################
chr_lengths <- c(
  "1" = 195471971, "2" = 182113224, "3" = 160039680, "4" = 156508116, 
  "5" = 151834684, "6" = 149736546, "7" = 145441459, "8" = 129401213, 
  "9" = 124595110, "10"= 130694993, "11"= 122082543, "12"= 120129022,
  "13"= 120421639, "14"= 124902244, "15"= 104043685, "16"= 98207768,  
  "17"= 94987271,  "18"= 90702639,  "19"= 61431566
)

gap_manhattan <- 7e6 
chr_offsets_manhattan <- cumsum(c(0, chr_lengths[-length(chr_lengths)] + gap_manhattan))
names(chr_offsets_manhattan) <- names(chr_lengths)
chr_mids_manhattan <- chr_offsets_manhattan + chr_lengths / 2

chr_bounds_manhattan <- tibble(
  Chr = names(chr_lengths),
  xmin = chr_offsets_manhattan,
  xmax = chr_offsets_manhattan + chr_lengths,
  shade = as.integer(names(chr_lengths)) %% 2 == 0
)

plot_df <- mpqtl_topLOD %>%
  mutate(
    Gene = as.character(Gene),
    Chr = as.character(Chr),
    abs_pos = Pos_bp + chr_offsets_manhattan[Chr],
    is_nash = tolower(Gene) %in% nash_list
  )

nash_genes <- plot_df %>%
  filter(is_nash) %>%
  distinct(Gene, Gene_chr, Gene_start) %>%
  mutate(
    Gene_chr = as.character(Gene_chr),
    gene_abs = Gene_start + chr_offsets_manhattan[Gene_chr]
  ) %>%
  arrange(Gene)

gene_colors <- setNames(scales::hue_pal()(nrow(nash_genes)), nash_genes$Gene)

ggplot() +
  geom_rect(data = chr_bounds_manhattan, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = shade), alpha = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "grey90", "FALSE" = "white")) +
  geom_vline(data = nash_genes, aes(xintercept = gene_abs, color = Gene), linewidth = 0.4, alpha = 0.6, show.legend = FALSE) +
  geom_point(data = filter(plot_df, !is_nash), aes(x = abs_pos, y = Lod), color = "grey50", size = 0.8, alpha = 0.5) +
  geom_point(data = filter(plot_df, is_nash), aes(x = abs_pos, y = Lod, color = Gene), size = 2, alpha = 1) +
  scale_color_manual(values = gene_colors, name = "NASH Gene") +
  scale_x_continuous(breaks = chr_mids_manhattan, labels = names(chr_mids_manhattan), expand = c(0.01, 0)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_text(size = 8), legend.key.height = unit(0.4, "cm"), legend.text = element_text(size = 7)) +
  labs(x = "Chromosome", y = "LOD Score", title = "Male pQTL LOD Scores")

################################################################################
# 6. Cis/Trans Density and Map Plots
################################################################################
gap_density <- 0 
chr_offsets_density <- cumsum(c(0, chr_lengths[-length(chr_lengths)] + gap_density))
names(chr_offsets_density) <- names(chr_lengths)
chr_mids_density <- chr_offsets_density + chr_lengths / 2

chr_bounds_density <- tibble(
  Chr = names(chr_lengths),
  xmin = chr_offsets_density,
  xmax = chr_offsets_density + chr_lengths,
  shade = as.integer(names(chr_lengths)) %% 2 == 0
)

plot_df_density <- mpqtl_topLOD %>%
  mutate(
    Chr = as.character(Chr),
    abs_pos = Pos_bp + chr_offsets_density[Chr],
    abs_gene = Gene_start + chr_offsets_density[as.character(Gene_chr)]
  )

x_scale <- scale_x_continuous(breaks = chr_mids_density, labels = names(chr_mids_density), expand = c(0.01, 0))

p_density <- ggplot() +
  geom_density(data = filter(plot_df_density, qtl_type == "trans"), aes(x = abs_pos), color = "firebrick", linewidth = 0.8, adjust = 0.25) +
  geom_density(data = filter(plot_df_density, qtl_type == "cis"), aes(x = abs_pos), color = "steelblue", linewidth = 0.8, adjust = 0.25) +
  x_scale +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(y = "Density")

p_map <- ggplot() +
  geom_rect(data = chr_bounds_density, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = shade), alpha = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "grey90", "FALSE" = "white")) +
  geom_point(data = plot_df_density, aes(x = abs_pos, y = abs_gene, color = qtl_type), size = 2, alpha = 0.7) +
  x_scale +
  scale_y_continuous(breaks = chr_mids_density, labels = names(chr_mids_density), expand = c(0.01, 0)) +
  scale_color_manual(values = c("cis" = "steelblue", "trans" = "firebrick"), name = "QTL Type") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) +
  labs(x = "QTL Position", y = "Gene Position")

p_density / p_map + plot_layout(heights = c(1, 4))