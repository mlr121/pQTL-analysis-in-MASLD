# pQTL-analysis-in-MASLD
 
Analysis code for proteome-wide QTL (pQTL) mapping in a genetically diverse MASLD mouse model (CC founders / F2 cross), with mediation analysis, ubiquitination (di-glycine) site profiling, and cross-species (mouse-to-human) validation of candidate proteins linking proteostasis to hepatic lipid metabolism. This code supports the manuscript *"Proteome-wide QTL mapping enables gene-protein-phenotype metabolic network construction in a genetically diverse MASLD mouse model."*

Many figures created within these scripts are also available in interactive format online at https://mlr98.shinyapps.io/pQTL_MASLD/
 
## Overview
 
The pipeline contains four sections:
 
1. **pQTL mapping** (01,02) — load LOD-score QTL results, annotate cis/trans status, and visualize the pQTL landscape (male and female mice analyzed in parallel, code is contained in separate files).
2. **Mediation analysis** (03) — test whether specific proteins (e.g., Ubr1) mediate the effect of a chromosome 2 QTL locus on hepatic steatosis.
3. **Ubiquitination site analysis** (04) — process di-glycine (GlyGly) enriched PTM data, normalize ubiquitination sites to protein abundance, and test for UBR1-dependent effects on steatosis.
4. **Human validation** (05) — process human liver proteomics/PTM data, map mouse-to-human orthologs, and test for concordance between mouse and human tau (association strength) statistics.

## Repo structure
```text
pQTL-analysis-in-MASLD/
├── README.md
├── 01_pqtl_proteomics_female.R
├── 02_pqtl_proteomics_male.R
├── 03_mediation_analysis.R
├── 04_ptm_ubiquitination_analysis.R
├── 05_human_cross_species_analysis.R
├── functions/
│   ├── gam_helpers.R
│   ├── mediation_helpers.R
│   └── impute_data.R
└── data/
    ├── raw/
      ├── mouse_liver_proteomics.zip/
       ├── Auwerx_Liver_4databasesearch_PG_FINAL
      ├── rawdata.zip/
       ├── final_Gene_Female_log2_filter05_ProMatrix_QTL.tsv
       ├── fixcol_ProMatrix_metadata.csv
       ├── F2_phenotypes.csv
       ├── mediation_analysis_marker.csv
       ├── pqtl_LIV_Stea.csv
       ├── simple_UbReport.txt
       ├── Mus_musculus_c57bl6nj.C57BL_6NJ_v1.pep.all.fasta
       ├── Mus_musculus_casteij.CAST_EiJ_v1.pep.all.fasta
       ├── Mus_musculus_pwkphj.PWK_PhJ_v1.pep.all.fasta
       ├── Mus_musculus_129s1svimj.129S1_SvImJ_v1.pep.all.fasta
      ├── human.zip/
           ├── human_proteomics.csv
           ├── human_metadata.csv
           └── 20260108_102414_HumanLiver_Report_PTM.csv
```
## Scripts
 
| Script | Purpose |
|---|---|
| `01_pqtl_proteomics_female.R` | Loads female mouse pQTL results, annotates cis/trans QTLs, generates Manhattan and cis/trans map plots, loads the female proteomics matrix and phenotype table, and computes protein–phenotype associations (GAM-based Kendall's tau) genome-wide and for Ubr1 specifically. |
| `02_pqtl_proteomics_male.R` | Same pipeline as above, run on the male mouse cohort, to enable sex-stratified comparison of pQTL LOD scores (male vs. female). |
| `03_mediation_analysis.R` | Runs causal mediation analysis (`mediation` package) testing whether individual proteins mediate the effect of chromosome 2 genotype markers on hepatic steatosis; includes a single-marker analysis fixed at the chr2 QTL peak, a genome-wide marker × protein sweep, and an Ubr1-as-mediator sweep across trans-associated proteins. |
| `04_ptm_ubiquitination_analysis.R` | Processes di-glycine (ubiquitination site) PTM data, normalizes site-level abundance to total protein abundance, tests Ubr1- and steatosis-associated ubiquitination sites, performs N-terminal/N-end rule classification from strain-specific proteome FASTAs, runs attenuation analysis (does adjusting for Ubr1 attenuate a site's effect on steatosis?), and runs GSEA/ORA pathway enrichment on ubiquitination tau statistics. |
| `05_human_cross_species_analysis.R` | Processes human liver proteomics and metadata, tests UBR1 abundance across diagnosis groups (e.g., NASH vs. not-NAFLD), computes protein–phenotype tau statistics, maps mouse-to-human orthologs via `biomaRt`, joins mouse and human tau statistics to test sign concordance, processes human di-glycine PTM data, and runs GSEA/ORA on human ubiquitination tau statistics. |
 
## Data files
 
These are the raw data files referenced by the scripts above.

**Mouse pQTL / proteomics**
- `final_Gene_Female_log2_filter05_ProMatrix_QTL.tsv` — female mouse pQTL mapping results (LOD scores, position, status)
- `fixcol_ProMatrix_metadata.csv` — mouse liver protein abundance matrix
- `F2_phenotypes.csv` — F2 cross mouse phenotype/trait table (includes Sex, Steatosis, etc.)
- `mediation_analysis_marker.csv` — genotype marker data used for mediation analysis (chr2 markers)
- `pqtl_LIV_Stea.csv` — liver steatosis phenotype table used in mediation analysis
- `simple_UbReport.txt` — mouse di-glycine (ubiquitination site) PTM report
  
**Mouse reference genomes (strain-specific proteomes, used for N-terminal/N-end rule analysis)**
- `Mus_musculus_c57bl6nj.C57BL_6NJ_v1.pep.all.fasta`
- `Mus_musculus_casteij.CAST_EiJ_v1.pep.all.fasta`
- `Mus_musculus_pwkphj.PWK_PhJ_v1.pep.all.fasta`
- `Mus_musculus_129s1svimj.129S1_SvImJ_v1.pep.all.fasta`
  
**Human liver data**
- `human_proteomics.csv` — human liver proteome (DIA report format)
- `human_metadata.csv` — human sample metadata (diagnosis, NAS score, medications, steatosis grade, etc.)
- `20260108_102414_HumanLiver_Report_PTM.csv` — human liver di-glycine PTM report
 
 
## Requirements
# Install CRAN dependencies
install.packages(c(
  "dplyr", "tidyr", "tibble", "readr", "data.table", "stringr",
  "ggplot2", "plotly", "patchwork", "ggsignif", "ggpmisc", "viridis",
  "scales", "broom", "purrr", "here", "mgcv", "mediation", "BioVenn"
))

# Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(c(
  "ensembldb", "EnsDb.Mmusculus.v79", "biomaRt",
  "msigdbr", "fgsea", "clusterProfiler", "Biostrings"
))

 
## Recommended Execution Order
Scripts depend on intermediate data objects created in previous stages (e.g., pqtl_topLOD and mpqtl_topLOD generated in scripts 01 and 02 are required for mediation and cross-species translation in scripts 03–05).

Execute the scripts sequentially in an active R session or workflow runner:

Female pQTL Analysis: 01_pqtl_proteomics_female.R

Male pQTL Analysis: 02_pqtl_proteomics_male.R

Causal Mediation Analysis: 03_mediation_analysis.R

Ubiquitination & N-End Rule Profiling: 04_ptm_ubiquitination_analysis.R

Human Translation & Ortholog Mapping: 05_human_cross_species_analysis.R

## Citation
 
If you use this code, please cite: *https://doi.org/10.64898/2026.07.26.740836*
 
## Contact
 
Meg Robinson, Coon Lab, University of Wisconsin–Madison mlrobinson3@wisc.edu
