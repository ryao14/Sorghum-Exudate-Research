# Sorghum-Exudate-Research

**Overview**
- Code accompanying my SRI research project on how Sorghum bicolor root exudates modulate the microbiome under high salinity stress.
- R scripts for sorghum rhizosphere/endosphere microbiome analysis: differential KEGG pathways, pathway ordination, and amplicon-based community summaries.
- Paths are hard-coded for Windows; update them to your locations before running.

**Requirements**
- R (>= 4.2) with: data.table, R.utils, phyloseq, ALDEx2, tidyverse, ggpicrust2, ggplot2, ggh4x, ggprism, vegan, tibble, edgeR, KEGGREST, stringr, psych, dplyr, readxl, reshape2, RColorBrewer, scales.
- BiocManager to install Bioconductor packages (phyloseq, KEGGREST, edgeR) if missing.
- Input files: PICRUSt2 KO and pathway outputs, sample metadata CSV, sequencing metadata TXT, and a Qiita-style Excel export.

**Data Inputs (edit paths in scripts)**
- PICRUSt2 KO metagenome table: picrust2_out/KO_metagenome_out/pred_metagenome_unstrat_cleaned.tsv.
- PICRUSt2 pathway abundance table: picrust2_out/pathways_out/path_abun_unstrat.tsv.
- Metadata: metadata.csv with sample_name and sample_type; sequencing metadata 15207_20250224-075556 (1).txt.
- Qiita export for OTUs and taxonomy: qiita_export1.xlsx.

**Scripts**
- Differential KEGG pathways: [edger_kegg_pathways.Rmd](edger_kegg_pathways.Rmd). Builds a KEGG abundance matrix from PICRUSt2 KO output, aligns with metadata, runs edgeR on endosphere vs rhizosphere samples, annotates pathways with ggpicrust2, and saves endosphere_vs_rhizosphere_plot.tiff.
- Pathway PCoA: [pathway_pcoa.Rmd](pathway_pcoa.Rmd). Loads PICRUSt2 pathway abundances, matches metadata, computes Bray-Curtis distances, performs PCoA, groups samples by the final digit of the ID (Plant 1-3), and saves pcoa_picrust2_pathways.png plus pcoa_results.csv.
- Amplicon phyloseq analyses: [phyloseq_genomic_analyses.r](phyloseq_genomic_analyses.r). Creates a phyloseq object from qiita_export1.xlsx, runs alpha diversity, Bray-Curtis PCoA with PERMANOVA, class/genus relative abundance summaries, Pearson correlations between taxa and soil metadata, core microbiome detection, clustering, and exports multiple CSV/TIFF outputs alongside OTU_phyloseq.rds.

**How to Run**
- Open the scripts in R or RStudio, update setwd and file paths, then run interactively. For Rmd files you can knit or render with rmarkdown::render("edger_kegg_pathways.Rmd") and rmarkdown::render("pathway_pcoa.Rmd").
- From a shell after setting paths: `Rscript phyloseq_genomic_analyses.r`.
- Confirm outputs in your working directory: PNG/TIFF pathway plots, alpha diversity tables, PERMANOVA results, relative abundance CSVs, correlation matrices, core taxa exports, and OTU_phyloseq.rds.

**Notes**
- Sample grouping assumes IDs end with .1/.2/.3 to denote Plant 1-3.
- Differential analysis filters to sample_type endosphere and rhizosphere; adjust the subset for other comparisons.
