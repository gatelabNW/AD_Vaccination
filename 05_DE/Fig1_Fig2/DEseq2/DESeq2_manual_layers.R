# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                         AN1792 Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 04-03-2024
# Written by: Anne Forsyth
# Summary: Differential Expression with DESeq2 for manual layers (Fig 1F-H, Ext Fig 1E, Ext Fig 2C)
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("DESeq2")
  library("ggrepel")
  library("UpSetR")
  library("randomcoloR")
  library("extrafont")
  library("showtext")
})

# Define output folder
output_folder <- "/path/to/manual/layer/deseq/results/"

# Load integrated Seurat object 
s <- readRDS("/path/to/all_samples_03.rds")

#-------------------------------------------------------------------------------
# Run differential expression

# Combine sample-level data
DefaultAssay(s) <- "Spatial"
s <- JoinLayers(s)

# Extract raw counts for each layer for each comparison, and store in a list
clusters <- unique(s@meta.data$manual_annotation)
comparison_counts <- list()
comparisons <- c("iAD_vs_nAD", "nAD_vs_NNC", "lim_vs_nAD", "ext_vs_nAD", "ext_vs_lim")
for (comparison in comparisons) {
  
  # Define idents for comparison
  ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
  ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
  
  # Generate a list of count matrices for each layer for the current comparison
  clust_counts <- list()
  for (cluster in clusters) {
    print(cluster)
    
    # Subset for current layer  
    clust_s <- subset(s, manual_annotation == cluster)
    
    # Extract raw counts
    counts <- GetAssayData(clust_s, assay = "Spatial", layer = "counts")
    
    # If comparing iAD-Extensive or iAD-Limited, set idents to condition_clearance
    if ((ident.1 %in% c("ext", "lim")) | (ident.2 %in% c("ext", "lim"))) {
      Idents(clust_s) <- "condition_clearance"
    } else {
      Idents(clust_s) <- "condition"
    }
    
    # Calculate percent expression for all genes based on raw counts 
    LFC <- FoldChange(object = clust_s, ident.1 = ident.1, ident.2 = ident.2, fc.name = "avg_log2FC", base = 2, layer = "counts") %>% data.frame()
    
    # Filter for genes expressed in at least 1% of either group
    genes_keep <- row.names(LFC)[LFC$pct.1 >= 0.01 | LFC$pct.2 >= 0.01]
    
    # Remove contamination
    if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
      genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)] 
    }
    
    # Filter count matrix
    clust_counts[[cluster]] <- counts[row.names(counts) %in% genes_keep,]
  }
  comparison_counts[[comparison]] <- clust_counts
}

# Extract metadata for each layer, and store in a list (meta data will be the same for all comparisons)
clust_meta <- list()
for (cluster in clusters) {
  clust_meta[[cluster]] <- s@meta.data[s@meta.data$manual_annotation == cluster,]
}

# Confirm raw count columns and metadata rows match for each comparison and each layer
for (comp in comparisons) {
  clust_counts <- comparison_counts[[comp]]
  for (cluster in clusters) {
    print(sum(colnames(clust_counts[[cluster]]) != row.names(clust_meta[[cluster]])))
  }
}

# Sum raw counts by sample ID to generate pseudobulk data
comparison_bulk_counts <- list()
for (comp in comparisons) {
  clust_counts <- comparison_counts[[comp]]
  bulk_counts <- list()
  for (cluster in clusters) {
    colnames(clust_counts[[cluster]]) <- clust_meta[[cluster]]$sample_id
    bulk <- t(rowsum(t(clust_counts[[cluster]]), group = colnames(clust_counts[[cluster]])))
    bulk_counts[[cluster]] <- bulk
  }
  comparison_bulk_counts[[comp]] <- bulk_counts
}

# Subset metadata to have one row per sample, and standardize continuous covariates 
for (cluster in clusters) {
  meta <- clust_meta[[cluster]][,c('age', 'sex', 'sample_id', 'gw_avg_features', 'gDNA_percent', 'condition', 'condition_clearance')]
  meta <- unique(meta)
  meta$age_centered <- as.numeric((meta$age - mean(meta$age)) / sd(meta$age))
  meta$gw_features_centered <- as.numeric((meta$gw_avg_features - mean(meta$gw_avg_features)) / sd(meta$gw_avg_features))
  meta$gDNA_centered <- as.numeric((meta$gDNA_percent - mean(meta$gDNA_percent)) / sd(meta$gDNA_percent))
  meta$sex <- as.factor(meta$sex)
  row.names(meta) <- meta$sample_id
  clust_meta[[cluster]] <- meta
}

# Confirm pseudobulk columns and metadata rows match for each comparison and each layer
for (comp in comparisons) {
  bulk_counts <- comparison_bulk_counts[[comp]]
  for (cluster in names(bulk_counts)) {
    meta <- clust_meta[[cluster]]
    bulk <- bulk_counts[[cluster]]
    bulk <- bulk[, row.names(meta)]
    
    print(sum(colnames(bulk) != row.names(meta)))
    print(sum(is.na(bulk)))
    
    bulk_counts[[cluster]] <- bulk
  }
  comparison_bulk_counts[[comp]] <- bulk_counts
}

# Run DESeq2 for each comparison, and save results in a list
comparison_results <- list()
for (comparison in comparisons) {
  bulk_counts <- comparison_bulk_counts[[comparison]]
  
  # Used to extract results later
  comp_name <- paste0("condition_", comparison)
  
  # Define idents for current comparison
  ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
  ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
  
  # Run DESeq2 for each layer
  clust_results <- list()
  for (cluster in clusters) {
    meta <- clust_meta[[cluster]]
    
    # If comparing iAD-Extensive or iAD-Limited, override condition variable
    if ((ident.1 %in% c("ext", "lim")) | (ident.2 %in% c("ext", "lim"))) {
      meta$condition <- factor(meta$condition_clearance)
    } else {
      meta$condition <- factor(meta$condition)
    }
    
    # Set condition reference level
    meta$condition <- relevel(meta$condition, ref = ident.2)
    
    # Variable used to skip to the next analysis if the current fails
    skip_to_next <- FALSE
    
    # Run DE using negative binomial GLM
    tryCatch(
      {
        dds <- DESeqDataSetFromMatrix(countData = bulk_counts[[cluster]], colData = meta, design= ~ sex + age_centered + gw_features_centered +
                                        gDNA_centered + condition)
        results <- DESeq(dds)
      },
      error = function(e) {
        skip_to_next <- TRUE
      }
    )
    if (skip_to_next) {
      print(paste0("Error in ", comparison, ", cluster ", cluster))
      next
    }
    
    tryCatch(
      {
        results <- results(results, name = paste0(comp_name))
        clust_results[[cluster]] <- results
      },
      error = function(e) {
        print(paste0("No ", comparison, " results for ", cluster))
      }
    )
  }
  comparison_results[[comparison]] <- clust_results
}

# Set DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Define DEGs 
for (comparison in comparisons) {
  clust_results <- comparison_results[[comparison]]
  for (cluster in names(clust_results)) {
    print(cluster)
    results <- clust_results[[cluster]]
    results$gene <- row.names(results)
    results$DE <- "Not DE"
    results$DE[results$log2FoldChange > fc_thresh & results$padj < p_thresh] <- "Upregulated"
    results$DE[results$log2FoldChange < -fc_thresh & results$padj < p_thresh] <- "Downregulated"
    
    # Create variable for volcano plot labels
    results$DE_gene <- NA
    results$DE_gene[results$DE != "Not DE"] <- results$gene[results$DE != "Not DE"]
    clust_results[[cluster]] <- results
  }
  comparison_results[[comparison]] <- clust_results
}

# Save results
for (comparison in comparisons) {
  clust_results <- comparison_results[[comparison]]
  for (cluster in names(clust_results)) {
    write.csv(clust_results[[cluster]], paste0(output_folder, comparison, "/deseq/results/", cluster, ".csv"))
  }
}

