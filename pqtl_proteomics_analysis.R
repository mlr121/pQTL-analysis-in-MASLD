#load data and filter to significant ----
fdata <- read.csv("rawdata/final_Gene_Female_log2_filter05_ProMatrix_QTL.tsv", sep = "\t", stringsAsFactors = F)


##filtering ----
library(dplyr)

fdata$Lod_threshold <- NULL

pqtls <- fdata %>%
  dplyr::filter(Status == "Significant_001" | Status == "Significant_005")

#annotate with cis and trans ----
library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(tidyr)

ensdb <- EnsDb.Mmusculus.v79
protein_locations <- genes(ensdb)
protein_locations_df <- as.data.frame(protein_locations)
protein_locations_df <- protein_locations_df[, c("symbol", "Gene_chr", "start", "end", "strand")]
protein_locations_df$Gene_chr <- as.factor(protein_locations_df$Gene_chr)
Chromosome_order <- c(as.character(1:19))  # This will create "1", "2", ..., "19"
Chromosome_order <- c(Chromosome_order, "X", "Y")
protein_locations_df$Gene_chr <- factor(protein_locations_df$Gene_chr, levels = Chromosome_order)
protein_locations_df$Gene_chr <- as.integer(protein_locations_df$Gene_chr)

pqtls$Phenotype <- gsub("_2", "", pqtls$Phenotype)
pqtls$Phenotype <- gsub("_3", "", pqtls$Phenotype)
pqtls$Phenotype <- gsub("_4", "", pqtls$Phenotype)
pqtls$Phenotype <- gsub("_5", "", pqtls$Phenotype)

pqtls_annotated <- pqtls %>%
  left_join(protein_locations_df, by = c("Phenotype" = "symbol"))

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

pqtl_topLOD <- pqtls_annotated[, .SD[which.max(Lod)], by = .(Phenotype, Chr, qtl_type)]
pqtl_topLOD<- pqtl_topLOD[order(-pqtl_topLOD$Lod), ]



qtl_counts <- pqtl_topLOD %>%
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
pqtl_topLOD <- pqtl_topLOD %>%
  mutate(mash_protein = as.factor(ifelse(tolower(Phenotype) %in% tolower(nash_list), Phenotype, NA)))


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

gap <- 7e6  # 3 Mb gap between chromosomes, adjust to taste

chr_offsets <- cumsum(c(0, chr_lengths[-length(chr_lengths)] + gap))
names(chr_offsets) <- names(chr_lengths)

chr_mids <- chr_offsets + chr_lengths / 2

chr_bounds <- tibble(
  Chr   = names(chr_lengths),
  xmin  = chr_offsets,
  xmax  = chr_offsets + chr_lengths,
  shade = as.integer(names(chr_lengths)) %% 2 == 0
)

plot_df <- pqtl_topLOD %>%
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
p <- ggplot(pqtl_topLOD, aes(x = qtl_loc, y = prot_loc, color = qtl_type, size = Lod, text = paste0(Phenotype, Chr))) +
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

plot_df <- pqtl_topLOD %>%
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
#compare male and female ----

male_selected <- mpqtl_topLOD[, c("Phenotype", "Chr", "qtl_type", "Lod")]
female_selected <- pqtl_topLOD[, c("Phenotype", "Chr", "qtl_type", "Lod", "nash_protein")]

mvf <- merge(male_selected, female_selected, by = c("Phenotype", "Chr", "qtl_type"), suffixes = c("_male", "_female"))

library(viridis)
library(ggpmisc)
mvf_scatter <- ggplot(mvf, aes(x = Lod_male, y = Lod_female, color = nash_protein)) +
  geom_point(aes(
    text = paste("Protein: ", Phenotype, "QTL Chromosome: ", Chr),
    size = ifelse(is.na(nash_protein), 0.5, 0.6),  # More subtle size difference
    alpha = ifelse(is.na(nash_protein), 2, 3)  # More subtle transparency difference
  )) +
  #scale_color_viridis(discrete = TRUE, na.value = "darkgrey") +  # Set NA values to grey
  labs(
    title = "Male v Female pQTL LOD scores",
    x = "Male",
    y = "Female"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(size = "none", alpha = "none") +  # Remove the size and alpha legends
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add best fit line (no standard error ribbon)
  stat_poly_eq(aes(label = paste("y = ", ..eq.label.., ", R² = ", ..rr.label..)), 
               label.x = 0.7, label.y = 0.8, color = "black", size = 5) 

ggplotly(mvf_scatter)

##venn diagram ----
library(BioVenn)

male_data <- apply(mpqtl_topLOD[, c("Phenotype", "Chr", "qtl_type")], 1, paste, collapse = "_")
female_data <- apply(pqtl_topLOD[, c("Phenotype", "Chr", "qtl_type")], 1, paste, collapse = "_")

listx <- male_data
listy <- female_data
listz <- NULL

draw.venn(listx, listy, listz, title = "shared QTLs", subtitle = NA, xtitle = "male", ytitle = "female")

#proteomics data load ----
proteomics <- read.csv(file="rawdata/fixcol_ProMatrix_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
row.names(proteomics) <- proteomics$Gene.Symbol

protein_matrix <- proteomics
protein_matrix$Gene.Symbol <- NULL
protein_matrix <- t(protein_matrix) 
protein_matrix <- as.data.frame(protein_matrix)%>%
  mutate(ID = row.names(protein_matrix))

pheno <- read.csv("rawdata/F2_phenotypes.csv", stringsAsFactors = F)
pheno <- pheno %>%
  mutate(ID = gsub("-", "_", ID))

non_traits <- c("ID","Strain","Generation","Sex")

pheno_cols <- colnames(pheno)[!(colnames(pheno) %in% non_traits)]
pheno_cols <- pheno_cols[sapply(pheno[pheno_cols], is.numeric)]

all_pheno <- protein_matrix %>%
  left_join(pheno, by = "ID")

#compute tau for all phenotypes ----
compute_gam_associations <- function(all_pheno, pheno_cols) {
  # Requires mgcv
  if (!"mgcv" %in% .packages(all.available = TRUE)) {
    # try to load; will error later if missing
    suppressPackageStartupMessages(require(mgcv))
  } else {
    suppressPackageStartupMessages(library(mgcv))
  }
  
  # Basic checks
  if (!all(pheno_cols %in% colnames(all_pheno))) {
    missing_cols <- pheno_cols[!pheno_cols %in% colnames(all_pheno)]
    stop("These pheno_cols are not in all_pheno: ", paste(missing_cols, collapse = ", "))
  }
  if (!"Sex" %in% colnames(all_pheno)) stop("all_pheno must contain a 'Sex' column.")
  
  # Ensure Sex is a factor (but keep original)
  sex_factor <- as.factor(all_pheno$Sex)
  
  # Identify numeric columns in the merged table
  numeric_cols <- colnames(all_pheno)[sapply(all_pheno, is.numeric)]
  
  # Proteins = numeric columns that are NOT the pheno_cols
  protein_cols <- setdiff(numeric_cols, pheno_cols)
  
  if (length(protein_cols) == 0) stop("No protein columns detected. Check pheno_cols selection.")
  
  results_list <- vector("list", length = length(pheno_cols))
  names(results_list) <- pheno_cols
  
  # Loop through phenotype columns
  for (trait in pheno_cols) {
    message("Processing trait: ", trait)
    y <- all_pheno[[trait]]
    
    # For each protein, compute GAM + tau
    trait_res <- lapply(protein_cols, function(prot) {
      protein_vals <- all_pheno[[prot]]
      
      # Require >= 3 paired non-missing observations
      if (sum(!is.na(protein_vals) & !is.na(y)) < 3) return(NULL)
      
      df_model <- data.frame(
        trait = y,
        sex = sex_factor,
        protein = protein_vals
      )
      
      fit <- tryCatch({
        mgcv::gam(trait ~ s(protein) + sex, data = df_model)
      }, error = function(e) return(NULL))
      
      if (is.null(fit)) return(NULL)
      
      s_summary <- summary(fit)
      
      # safe extraction of p-value (in case s.table missing)
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
      
      direction <- if (is.na(tau_val)) NA_character_ else if (tau_val > 0) "1" else if (tau_val < 0) "-1" else "neutral"
      
      data.frame(
        Protein = prot,
        p_value = pval,
        deviance_explained = dev_expl,
        tau = tau_val,
        direction = direction,
        stringsAsFactors = FALSE
      )
    })
    
    # bind and remove NULLs (returned when <3 obs or fit failed)
    trait_df <- do.call(rbind, trait_res)
    if (!is.null(trait_df) && nrow(trait_df) > 0) {
      results_list[[trait]] <- trait_df
    } else {
      results_list[[trait]] <- NULL
    }
  }
  
  return(results_list)
}

##run tau calc ----
results_tau_allpheno <- compute_gam_associations(all_pheno, pheno_cols)

##combine ----
combined_df <- do.call(rbind, lapply(names(results_tau_allpheno), function(trait) {
  df <- results_tau_allpheno[[trait]]
  if (!is.null(df)) df$Trait <- trait
  df
}))

#tau for ubr1 ----
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

ubr1_results <- compute_ubr1_associations(all_pheno, ubr1_col = "Ubr1")


