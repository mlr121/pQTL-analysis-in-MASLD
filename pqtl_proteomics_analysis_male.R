#load data and filter to significant ----
mdata <- read.csv("rawdata/final_Gene_Male_log2_filter05_ProMatrix_QTL.tsv", sep = "\t", stringsAsFactors = F)

##filtering ----
library(dplyr)

mdata$Lod_threshold <- NULL

pqtls <- mdata %>%
  dplyr::filter(Status == "Significant_001" | Status == "Significant_005")

#annotate with cis and trans ----
library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(tidyr)

##annotating ----
cis_distance_threshold <- 1000000

pqtls_annotated <- pqtls %>%
  mutate(
    qtl_type = case_when(
      # Cis: same Chromosome and within the distance threshold
      (Gene_chr == Chr) & (abs(Gene_start - Pos_bp) <= cis_distance_threshold) ~ "cis",
      
      # Trans: different Chromosomes or outside the distance threshold
      (Gene_chr != Chr) | (abs(Gene_start - Pos_bp) > cis_distance_threshold) ~ "trans"
    )
  )

#filter to only top LOD point for each pQTL ----
library(data.table)
library(ggplot2)
library(plotly)
setDT(pqtls_annotated)

mpqtl_topLOD <- pqtls_annotated[, .SD[which.max(Lod)], by = .(Phenotype, Chr, qtl_type)]
mpqtl_topLOD<- mpqtl_topLOD[order(-mpqtl_topLOD$Lod), ]



qtl_counts <- mpqtl_topLOD %>%
  count(qtl_type)

print(qtl_counts)

#annotate proteins of interest ----
nash_list <- c("Pnpla3", "Tm6sf2", "Gckr", "Mboat7", "Hsd17b13__1", "Fasn", "Acaca", "Acly",
               "Scd1", "Srebf1", "Mlxipl", "Hmgcs1", "Hmgcr", "Dgat2", "Nr1h4", "Cd36",
               "Fabp1", "Fabp4", "Lpl", "Ldlr", "Plin2", "Apoe", "Apob", "Cpt1a", "Acox1",
               "Cyp2e1", "Ppara", "Pparg", "Thrb", "Hnf4a", "Pck1", "Col1a1", "Acta2",
               "Tgfb1", "Timp1", "Mmp2", "Mmp9", "Tnf", "Trem2", "Socs3", "Il10",
               "Fgf21", "Gdf15", "Adipoq", "Saa1", "Sod2", "Gpx1", "Hmox1", "Gpam",
               "Pnpla2"
)

mpqtl_topLOD <- mpqtl_topLOD %>%
  mutate(mash_protein = as.factor(ifelse(
    sapply(tolower(Phenotype), function(p) any(str_detect(p, fixed(tolower(nash_list))))),
    Phenotype, 
    NA
  )))

#plot manhattan ----
# Mouse chr sizes (mm10/GRCm38)
chr_lengths <- c(
  "1"  = 195471971, "2"  = 182113224, "3"  = 160039680,
  "4"  = 156508116, "5"  = 151834684, "6"  = 149736546,
  "7"  = 145441459, "8"  = 129401213, "9"  = 124595110,
  "10" = 130694993, "11" = 122082543, "12" = 120129022,
  "13" = 120421639, "14" = 124902244, "15" = 104043685,
  "16" = 98207768,  "17" = 94987271,  "18" = 90702639,
  "19" = 61431566
)

gap <- 7e6 

chr_offsets <- cumsum(c(0, chr_lengths[-length(chr_lengths)] + gap))
names(chr_offsets) <- names(chr_lengths)

chr_mids <- chr_offsets + chr_lengths / 2

chr_bounds <- tibble(
  Chr   = names(chr_lengths),
  xmin  = chr_offsets,
  xmax  = chr_offsets + chr_lengths,
  shade = as.integer(names(chr_lengths)) %% 2 == 0
)

plot_df <- mpqtl_topLOD %>%
  mutate(
    Gene    = as.character(Gene),
    Chr     = as.character(Chr),
    abs_pos = Pos_bp + chr_offsets[Chr],
    is_nash = Gene %in% nash_list
  )

nash_genes <- plot_df %>%
  filter(is_nash) %>%
  distinct(Gene, Gene_chr, Gene_start) %>%
  mutate(
    Gene_chr = as.character(Gene_chr),
    gene_abs = Gene_start + chr_offsets[Gene_chr]
  ) %>%
  arrange(Gene)

gene_colors <- setNames(
  scales::hue_pal()(nrow(nash_genes)),
  nash_genes$Gene
)
ggplot() +
  
  geom_rect(
    data = chr_bounds,
    aes(xmin = xmin, xmax = xmax,
        ymin = -Inf, ymax = Inf,
        fill = shade),
    alpha = 0.3,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = c("TRUE" = "grey90", "FALSE" = "white")
  ) +
  
  geom_vline(
    data = nash_genes,
    aes(xintercept = gene_abs, color = Gene),
    linewidth = 0.4,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  
  geom_point(
    data = filter(plot_df, !is_nash),
    aes(x = abs_pos, y = Lod),
    color = "grey50",
    size = 0.8,
    alpha = 0.5
  ) +
  
  geom_point(
    data = filter(plot_df, is_nash),
    aes(x = abs_pos, y = Lod, color = Gene),
    size = 2,
    alpha = 1
  ) +
  
  scale_color_manual(
    values = gene_colors,
    name = "NASH Gene"
  ) +
  
  scale_x_continuous(
    breaks = chr_mids,
    labels = names(chr_mids),
    expand = c(0.01, 0)
  ) +
  
  theme_minimal() +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8),
    legend.key.height = unit(0.4, "cm"),
    legend.text = element_text(size = 7)
  ) +
  
  labs(
    x = "Chromosome",
    y = "LOD Score",
    title = "pQTL LOD Scores"
  )

#plot cis and trans ----
p <- ggplot(mpqtl_topLOD, aes(x = qtl_loc, y = prot_loc, color = qtl_type, size = Lod, text = paste0(Phenotype, Chr))) +
  geom_point() +
  scale_size_continuous(range = c(0.1,5)) +
  scale_color_manual(values = c("red", "blue")) +
  scale_x_continuous(breaks = seq(0, 2e9, by = 1e8)) +
  labs(
    title = "Protein Gene Locations vs. pQTL Locations",
    x = "Chromosome",
    y = "Protein Gene Location") + 
  theme_minimal() 

ggplotly(p)

gap <- 0  # 0 Mb gap between chromosomes, adjust to taste

chr_offsets <- cumsum(c(0, chr_lengths[-length(chr_lengths)] + gap))
names(chr_offsets) <- names(chr_lengths)

chr_mids <- chr_offsets + chr_lengths / 2

chr_bounds <- tibble(
  Chr   = names(chr_lengths),
  xmin  = chr_offsets,
  xmax  = chr_offsets + chr_lengths,
  shade = as.integer(names(chr_lengths)) %% 2 == 0
)

plot_df <- mpqtl_topLOD %>%
  mutate(
    Chr       = as.character(Chr),
    abs_pos   = Pos_bp + chr_offsets[Chr],
    abs_gene  = Gene_start + chr_offsets[as.character(Gene_chr)]
  )

library(patchwork)

# --- Shared x scale ---
x_scale <- scale_x_continuous(
  breaks = chr_mids,
  labels = names(chr_mids),
  expand = c(0.01, 0)
)

# --- Density plot ---
p_density <- ggplot() +
  geom_density(
    data = filter(plot_df, qtl_type == "trans"),
    aes(x = abs_pos),
    color = "firebrick", size = 0.8, adjust = 0.25
  ) +
  geom_density(
    data = filter(plot_df, qtl_type == "cis"),
    aes(x = abs_pos),
    color = "steelblue", size = 0.8, adjust = 0.25
  ) +
  x_scale +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(y = "Density")

# --- pQTL map ---
p_map <- ggplot() +
  geom_rect(
    data = chr_bounds,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = shade),
    alpha = 0.3, show.legend = FALSE
  ) +
  scale_fill_manual(values = c("TRUE" = "grey90", "FALSE" = "white")) +
  geom_point(
    data = plot_df,
    aes(x = abs_pos, y = abs_gene, color = qtl_type),
    size = 2, alpha = 0.7
  ) +
  x_scale +
  scale_y_continuous(
    breaks = chr_mids,
    labels = names(chr_mids),
    expand = c(0.01, 0)
  ) +
  scale_color_manual(
    values = c("cis" = "steelblue", "trans" = "firebrick"),
    name = "QTL Type"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(x = "QTL Position", y = "Gene Position")

# --- Combine ---
p_density / p_map + plot_layout(heights = c(1, 4))
