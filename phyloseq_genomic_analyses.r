library(phyloseq)
library(vegan)
library(ggplot2)
library(psych)
library(dplyr)
library(readxl)
library(reshape2)
library(RColorBrewer)

# Set working directory
setwd("C:/Users/ryany/Downloads/SRI study/study_15207_061225-081635/study_15207_061225-081635/BIOM/180200")

# Define column names
otu_col <- "otu"
sample_col <- "sample"

# Function to create phyloseq object from Excel file
create_phyloseq <- function(excel_file, otu_col = "otu", sample_col = "sample") {
  tryCatch({
    otu_df <- read_excel(excel_file, sheet = "OTU matrix")
    taxonomy_df <- read_excel(excel_file, sheet = "Taxonomy table")
    sample_df <- read_excel(excel_file, sheet = "Samples")
    
    otu_matrix <- as.matrix(otu_df[, setdiff(colnames(otu_df), otu_col)])
    rownames(otu_matrix) <- otu_df[[otu_col]]
    otu <- otu_table(otu_matrix, taxa_are_rows = TRUE)
    
    if (!all(c("otu", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") %in% colnames(taxonomy_df))) {
      stop("Taxonomy table must contain columns: otu, Kingdom, Phylum, Class, Order, Family, Genus, Species")
    }
    otu_ids <- rownames(otu_matrix)
    tax_ids <- taxonomy_df$otu
    common_ids <- intersect(otu_ids, tax_ids)
    if (length(common_ids) == 0) {
      stop("No matching OTU IDs between OTU matrix and taxonomy table.")
    }
    if (length(common_ids) < length(otu_ids)) {
      otu_matrix <- otu_matrix[common_ids, , drop = FALSE]
      otu <- otu_table(otu_matrix, taxa_are_rows = TRUE)
    }
    taxonomy_df <- taxonomy_df[taxonomy_df$otu %in% common_ids, ]
    invalid_phyla <- is.na(taxonomy_df$Phylum) | taxonomy_df$Phylum == "" | taxonomy_df$Phylum == "unclassified"
    invalid_genera <- is.na(taxonomy_df$Genus) | taxonomy_df$Genus == "" | taxonomy_df$Genus == "unclassified"
    if (sum(invalid_phyla) > 0) {
      taxonomy_df <- taxonomy_df[!invalid_phyla, ]
      otu_matrix <- otu_matrix[taxonomy_df$otu, , drop = FALSE]
      otu <- otu_table(otu_matrix, taxa_are_rows = TRUE)
    }
    if (sum(invalid_genera) > 0) {
      taxonomy_df <- taxonomy_df[!invalid_genera, ]
      otu_matrix <- otu_matrix[taxonomy_df$otu, , drop = FALSE]
      otu <- otu_table(otu_matrix, taxa_are_rows = TRUE)
    }
    tax_matrix <- as.matrix(taxonomy_df[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
    rownames(tax_matrix) <- taxonomy_df$otu
    tax <- tax_table(tax_matrix)
    
    if (!(sample_col %in% colnames(sample_df))) {
      stop("Sample column '", sample_col, "' not found in Samples sheet.")
    }
    sample_df <- data.frame(sample_df, check.names = FALSE)
    rownames(sample_df) <- sample_df[[sample_col]]
    sample <- sample_data(sample_df)
    
    sample_names_otu <- colnames(otu_matrix)
    sample_names_meta <- rownames(sample_df)
    if (!all(sample_names_otu %in% sample_names_meta)) {
      stop("Sample IDs in OTU matrix do not match those in Samples sheet.")
    }
    
    physeq <- phyloseq(otu, tax, sample)
    return(physeq)
  }, error = function(e) {
    return(NULL)
  })
}

# Excel file
excel_file <- "qiita_export1.xlsx"

# Create phyloseq object
physeq <- create_phyloseq(excel_file, otu_col = "otu", sample_col = "sample")
if (is.null(physeq)) {
  stop("Failed to create phyloseq object. Check Excel file structure and metadata.")
}

# Diversity Analyses
analyze_diversity <- function(physeq, name = "OTU") {
  tryCatch({
    if (dim(sample_data(physeq))[1] == 0) {
      stop("Sample data has zero dimensions. Check metadata in Samples sheet.")
    }
    
    alpha <- estimate_richness(physeq, measures = c("Shannon", "Observed"))
    alpha$SampleID <- rownames(alpha)
    write.csv(alpha, paste0(name, "_alpha_diversity.csv"), row.names = FALSE)
    
    dist <- phyloseq::distance(physeq, method = "bray")
    pcoa <- cmdscale(dist, k = 2, eig = TRUE)
    pcoa_df <- data.frame(SampleID = sample_names(physeq), PC1 = pcoa$points[,1], PC2 = pcoa$points[,2], 
                          Eigenvalue_PC1 = pcoa$eig[1], Eigenvalue_PC2 = pcoa$eig[2])
    sample_df <- as(sample_data(physeq), "data.frame")
    sample_df$SampleID <- rownames(sample_df)
    
    pcoa_df <- merge(pcoa_df, sample_df, by = "SampleID", all.x = TRUE)
    pcoa_df$PlantGroup <- substr(pcoa_df$SampleID, nchar(pcoa_df$SampleID), nchar(pcoa_df$SampleID))
    pcoa_df$Shape <- ifelse(pcoa_df$PlantGroup == "1", 15, ifelse(pcoa_df$PlantGroup == "2", 16, 17))
    
    # Calculate percentage of variance explained
    total_variance <- sum(pcoa$eig)
    percent_var <- 100 * pcoa$eig / total_variance
    pc1_percent <- round(percent_var[1], 1)
    pc2_percent <- round(percent_var[2], 1)
    
    p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = sample_type, shape = factor(PlantGroup))) +
      geom_jitter(size = 3, width = 0.05, height = 0.05) +
      scale_shape_manual(values = c("1" = 15, "2" = 16, "3" = 17), labels = c(".1" = ".1", ".2" = ".2", ".3" = ".3")) +
      theme_minimal() +
      labs(x = paste0("PC1 (", pc1_percent, "%)"), y = paste0("PC2 (", pc2_percent, "%)"), shape = "Plant Group") +
      guides(color = guide_legend("Sample Type"), shape = guide_legend("Plant Group")) +
      theme(
        text = element_text(family = "Times New Roman", color = "black"),
        axis.title = element_text(size = 24, color = "black"),  # Increased from 10
        axis.text = element_text(size = 18, color = "black"),   # Increased from 8
        legend.title = element_text(size = 24, color = "black"), # Increased from 10
        legend.text = element_text(size = 18, color = "black")   # Increased from 8
      )
    ggsave(paste0(name, "_pcoa_plot.tiff"), p, dpi = 300, device = "tiff", width = 8, height = 6, units = "in", compression = "lzw")
    
    permanova <- adonis2(dist ~ sample_type, data = sample_df)
    write.csv(as.data.frame(permanova), paste0(name, "_permanova.csv"), row.names = TRUE)
    
    return(list(alpha = alpha, pcoa = pcoa_df, permanova = permanova, eigenvalues = pcoa$eig))
  }, error = function(e) {
    return(NULL)
  })
}

diversity_results <- analyze_diversity(physeq)

# Taxonomic Profiling
summarize_taxonomy <- function(physeq, name = "OTU", level = "Class") {
  tryCatch({
    physeq_tax <- tax_glom(physeq, taxrank = level)
    tax_table <- as.data.frame(tax_table(physeq_tax))
    otu_table <- otu_table(physeq_tax, taxa_are_rows = TRUE)
    otu_table_df <- as.data.frame(otu_table)
    
    tax_names <- make.unique(as.character(tax_table[[level]]))
    rownames(otu_table_df) <- tax_names
    
    post_glom_counts <- data.frame(Taxon = tax_names, Sample_Count = rowSums(otu_table_df > 0))
    write.csv(post_glom_counts, paste0(name, "_post_glom_", level, "_summary.csv"), row.names = FALSE)
    
    rel_abund <- t(t(otu_table_df) / colSums(otu_table_df)) * 100
    write.csv(rel_abund, paste0(name, "_relative_abundance_", level, ".csv"), row.names = TRUE)
    
    top_taxa <- names(sort(rowSums(rel_abund), decreasing = TRUE))[1:min(10, nrow(rel_abund))]
    rel_abund_top <- rel_abund[top_taxa, , drop = FALSE]
    rel_abund_df <- reshape2::melt(as.matrix(rel_abund_top), varnames = c("Taxon", "SampleID"), value.name = "Abundance")
    
    sample_df <- as(sample_data(physeq), "data.frame")
    sample_df$SampleID <- rownames(sample_df)
    
    rel_abund_df <- merge(rel_abund_df, sample_df, by = "SampleID", all.x = TRUE)
    
    rel_abund_df <- rel_abund_df %>%
      mutate(SampleID = as.character(SampleID)) %>%
      mutate(SampleOrder = substr(SampleID, nchar(SampleID), nchar(SampleID))) %>%
      arrange(SampleOrder)
    
    p <- ggplot(rel_abund_df, aes(x = factor(SampleID, levels = unique(rel_abund_df$SampleID)), y = Abundance, fill = Taxon)) +
      geom_bar(stat = "identity", position = "stack") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Sample", y = "Relative Abundance (%)") +
      scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
      theme(text = element_text(family = "Times New Roman"),  # Set font family
            axis.title = element_text(size = 10),           # Axis labels at 10 pt
            axis.text = element_text(size = 8),             # Axis text at 8 pt
            legend.title = element_text(size = 10),         # Legend title at 10 pt
            legend.text = element_text(size = 8))           # Legend text at 8 pt
    ggsave(paste0(name, "_taxa_barplot_", level, ".tiff"), p, dpi = 300, device = "tiff", width = 12, height = 8, units = "in", compression = "lzw")
    
    return(rel_abund)
  }, error = function(e) {
    return(NULL)
  })
}

taxa_results_class <- summarize_taxonomy(physeq, "OTU", "Class")
taxa_results_genus <- summarize_taxonomy(physeq, "OTU", "Genus")

# Correlation with Metadata
compute_correlation <- function(physeq, name = "OTU", metadata_vars = c("soil_ph_1v1", "soi_h2o_total_organic_c", "soil_available_n", "ec_uspercm")) {
  tryCatch({
    physeq_class <- tax_glom(physeq, taxrank = "Class")
    otu_table <- as.data.frame(otu_table(physeq_class, taxa_are_rows = TRUE))
    metadata <- as.data.frame(sample_data(physeq))
    
    if (ncol(otu_table) != nrow(metadata)) {
      stop("Incompatible dimensions: OTU table has ", ncol(otu_table), " samples, but metadata has ", nrow(metadata), " samples.")
    }
    
    metadata <- metadata[, metadata_vars[metadata_vars %in% colnames(metadata)], drop = FALSE]
    metadata <- as.data.frame(lapply(metadata, as.numeric))
    metadata <- metadata[, colSums(is.na(metadata)) < nrow(metadata), drop = FALSE]
    
    if (ncol(metadata) == 0) {
      return(NULL)
    }
    
    rel_abund <- t(t(otu_table) / colSums(otu_table)) * 100
    mean_abund <- rowMeans(rel_abund)
    significant_classes <- names(mean_abund[mean_abund >= 1])
    
    if (length(significant_classes) > 0) {
      otu_table_signif <- otu_table[significant_classes, , drop = FALSE]
      other_abund <- colSums(otu_table[!(rownames(otu_table) %in% significant_classes), , drop = FALSE])
      otu_table <- rbind(otu_table_signif, other = other_abund)
    } else {
      other_abund <- colSums(otu_table)
      otu_table <- data.frame(t(other_abund), row.names = "other")
    }
    
    otu_table <- t(otu_table)
    corr_matrix <- corr.test(otu_table, metadata, method = "pearson")$r
    corr_pval <- corr.test(otu_table, metadata, method = "pearson")$p
    
    corr_df <- reshape2::melt(corr_matrix)
    pval_df <- reshape2::melt(corr_pval)
    corr_df$pval <- pval_df$value
    corr_df$Signif <- ifelse(corr_df$pval < 0.05, "*", "")
    
    tax_table <- as.data.frame(tax_table(physeq_class))
    corr_df$Var1 <- ifelse(corr_df$Var1 == "other", "other", tax_table[as.character(corr_df$Var1), "Class"])
    corr_df$Var1 <- as.factor(corr_df$Var1)
    
    p <- ggplot(corr_df, aes(x = Var2, y = Var1, fill = value, label = Signif)) +
      geom_tile() +
      geom_text() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
      theme_minimal() +
      labs(x = "Metadata", y = "Classes") +
      theme(axis.text.y = element_text(size = 6)) +  # Adjusted separately due to small text size
      theme(text = element_text(family = "Times New Roman"),  # Set font family
            axis.title = element_text(size = 10),           # Axis labels at 10 pt
            axis.text = element_text(size = 8),             # Axis text at 8 pt
            legend.title = element_text(size = 10),         # Legend title at 10 pt
            legend.text = element_text(size = 8))           # Legend text at 8 pt
    ggsave(paste0(name, "_correlation_heatmap.tiff"), p, dpi = 300, device = "tiff", width = 10, height = 8, units = "in", compression = "lzw")
    
    write.csv(corr_matrix, paste0(name, "_correlation_matrix.csv"), row.names = TRUE)
    return(corr_matrix)
  }, error = function(e) {
    return(NULL)
  })
}

corr_results <- compute_correlation(physeq)

# Core Microbiome Analysis
find_core_taxa <- function(physeq, name = "OTU", detection = 0.0001, prevalence = 0.1, abund_threshold = 1) {
  tryCatch({
    otu_table <- otu_table(physeq, taxa_are_rows = TRUE)
    if (any(dim(otu_table) == 0)) {
      stop("OTU table has zero dimensions.")
    }
    
    col_sums <- colSums(otu_table)
    sparsity_stats <- data.frame(Sample = colnames(otu_table), TotalCounts = col_sums, NonZeroOTUs = colSums(otu_table > 0))
    write.csv(sparsity_stats, paste0(name, "_sparsity_stats.csv"), row.names = FALSE)
    
    if (any(col_sums == 0)) {
      message("Warning: ", sum(col_sums == 0), " samples have zero total abundance.")
    }
    
    rel_abund <- t(t(otu_table) / pmax(col_sums, 1)) * 100
    prevalence <- apply(rel_abund >= detection, 1, mean)
    write.csv(data.frame(Taxon = names(prevalence), Prevalence = prevalence), paste0(name, "_prevalence_summary.csv"), row.names = FALSE)
    
    core_taxa_idx <- prevalence >= prevalence
    if (sum(core_taxa_idx) == 0) {
      core_taxa_idx <- prevalence >= 0.05
    }
    if (sum(core_taxa_idx) == 0) {
      stop("No taxa meet the relaxed prevalence threshold of 0.05.")
    }
    
    core_taxa_names <- names(which(core_taxa_idx))
    core_taxa <- prune_taxa(core_taxa_names, physeq)
    if (is.null(core_taxa)) {
      stop("Failed to prune taxa.")
    }
    
    core_taxa_class <- tax_glom(core_taxa, taxrank = "Class")
    core_abund <- as.data.frame(otu_table(core_taxa_class))
    core_tax <- as.data.frame(tax_table(core_taxa_class))
    write.csv(cbind(core_tax, core_abund), paste0(name, "_core_taxa_class.csv"), row.names = TRUE)
    
    core_rel_abund <- t(t(core_abund) / colSums(core_abund)) * 100
    mean_abund <- rowMeans(core_rel_abund)
    significant_classes <- rownames(core_tax)[mean_abund >= abund_threshold]
    other_abund <- colSums(core_rel_abund[!(rownames(core_rel_abund) %in% significant_classes), , drop = FALSE])
    
    core_df <- reshape2::melt(as.matrix(core_rel_abund[significant_classes, , drop = FALSE]), varnames = c("Taxon", "SampleID"), value.name = "Abundance")
    core_df$Taxon <- core_tax[match(core_df$Taxon, rownames(core_tax)), "Class"]
    core_df$Taxon <- ifelse(is.na(core_df$Taxon) | core_df$Taxon == "" | core_df$Taxon == "unclassified", "unknown", core_df$Taxon)
    # Remove "c_" and brackets from taxon names
    # Clean taxon names: remove "c_", square brackets, and any leading/trailing underscores or whitespace
    core_df$Taxon <- gsub("^c_", "", core_df$Taxon)                   # remove leading 'c_'
    core_df$Taxon <- gsub("\\[|\\]", "", core_df$Taxon)              # remove [ and ]
    core_df$Taxon <- gsub("^_+|_+$", "", core_df$Taxon)              # remove leading/trailing underscores
    core_df$Taxon <- gsub("^c__", "", core_df$Taxon)                 # in case there's a double underscore format
    core_df$Taxon <- trimws(core_df$Taxon)                           # remove leading/trailing whitespace
    
    
    other_df <- data.frame(Taxon = "other", SampleID = colnames(core_rel_abund), Abundance = other_abund)
    core_df <- rbind(core_df, other_df)
    
    if (any(is.na(core_df$Abundance))) {
      core_df <- core_df[!is.na(core_df$Abundance), ]
    }
    if (all(core_df$Abundance == 0)) {
      core_df$Abundance <- core_df$Abundance + 1e-10
    }
    
    sample_df <- as(sample_data(physeq), "data.frame")
    sample_df$SampleID <- rownames(sample_df)
    
    core_df <- merge(core_df, sample_df, by = "SampleID", all.x = TRUE)
    
    core_df <- core_df %>%
      mutate(SampleID = as.character(SampleID)) %>%
      mutate(SampleOrder = substr(SampleID, nchar(SampleID), nchar(SampleID))) %>%
      arrange(SampleOrder)
    
    # Define a custom order of distinct colors
    n_taxa <- length(unique(core_df$Taxon))
    base_colors <- brewer.pal(min(12, n_taxa), "Set3")
    custom_color_order <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", 
                            "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
    if (n_taxa > 12) {
      additional_colors <- rainbow(n_taxa - 12, s = 0.8, v = 0.9)
      all_colors <- c(custom_color_order[1:min(12, n_taxa)], additional_colors)
    } else {
      all_colors <- custom_color_order[1:n_taxa]
    }
    # Match colors to original taxon order
    unique_taxa <- unique(core_df$Taxon)
    all_colors <- all_colors[1:length(unique_taxa)]
    names(all_colors) <- unique_taxa
    
    # Define manual labels for 12 samples (4 compartments × 3 PlantGroups)
    manual_labels <- c("Bulksoil 1", "Endosphere 1", "Rhizoplane 1", "Rhizosphere 1",
                       "Bulksoil 2", "Endosphere 2", "Rhizoplane 2", "Rhizosphere 2",
                       "Bulksoil 3", "Endosphere 3", "Rhizoplane 3", "Rhizosphere 3")
    
    # Define manual labels for 12 samples (4 compartments × 3 PlantGroups)
    manual_labels <- c("Bulksoil 1", "Endosphere 1", "Rhizoplane 1", "Rhizosphere 1",
                       "Bulksoil 2", "Endosphere 2", "Rhizoplane 2", "Rhizosphere 2",
                       "Bulksoil 3", "Endosphere 3", "Rhizoplane 3", "Rhizosphere 3")
    
    p <- ggplot(core_df, aes(x = factor(SampleID, levels = unique(core_df$SampleID), labels = manual_labels), y = Abundance, fill = Taxon)) +
      geom_bar(stat = "identity", position = "stack") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right", legend.key.size = unit(0.5, "cm")) +
      labs(x = "Sample", y = "Relative Abundance (%)", fill = "Class") +
      scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
      scale_fill_manual(values = all_colors) +
      guides(fill = guide_legend(title = "Class")) +
      theme(
        text = element_text(family = "Times New Roman", color = "black"),
        axis.title = element_text(size = 10),  # Reverted to original 10 pt
        axis.text = element_text(size = 8),    # Reverted to original 8 pt
        legend.title = element_text(size = 10), # Reverted to original 10 pt
        legend.text = element_text(size = 8)    # Reverted to original 8 pt
      )
    ggsave(paste0(name, "_core_taxa_barplot_paper.tiff"), p, dpi = 300, device = "tiff", width = 12, height = 8, units = "in", compression = "lzw")
    core_rel_abund <- t(t(core_abund) / colSums(core_abund)) * 100
    dist_matrix <- vegdist(t(core_rel_abund), method = "bray")
    hclust_result <- hclust(dist_matrix, method = "ward.D2")
    cluster_groups <- cutree(hclust_result, k = 3)
    cluster_df <- data.frame(SampleID = names(cluster_groups), Cluster = as.factor(cluster_groups))
    
    sample_df <- as(sample_data(physeq), "data.frame")
    sample_df$SampleID <- rownames(sample_df)
    cluster_df <- merge(cluster_df, sample_df, by = "SampleID", all.x = TRUE)
    write.csv(cluster_df, paste0(name, "_sample_clusters.csv"), row.names = FALSE)
    
    return(core_taxa_class)
  }, error = function(e) {
    return(NULL)
  })
}

core_results <- find_core_taxa(physeq)

# Save phyloseq object
saveRDS(physeq, "OTU_phyloseq.rds")

cat("Analysis complete.\n")