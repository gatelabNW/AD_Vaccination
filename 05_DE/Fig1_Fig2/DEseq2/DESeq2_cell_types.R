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
# Summary: Differential Expression with DESeq2 for C2L enriched cell types (Fig 1M-N, Ext Fig 1F-G)
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
})

# Define output folder
output_folder <- "/path/to/celltype/deseq/results/"

# Load integrated Seurat object 
s <- readRDS("/path/to/all_samples_03.rds")

# Define filter operator 
`%notin%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
# Run differential expression

# Combine sample-level data
DefaultAssay(s) <- "Spatial"
s <- JoinLayers(s)

# Extract raw counts for each cell type for each comparison, and store in a list
cell_types <- colnames(s@meta.data)[34:51]
comparison_counts <- list()
comparisons <- c("iAD_vs_nAD", "nAD_vs_NNC", "lim_vs_nAD", "ext_vs_nAD", "ext_vs_lim")
for (comparison in comparisons) {
  
  # Define idents for comparison
  ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
  ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
  
  # Generate a list of count matrices for each cell type for the current comparison
  type_counts <- list()
  for (type in cell_types) {
    print(type)
    
    # Set Idents to cell type column to allow filtering using subset function 
    Idents(s) <- type
    
    # If there are no spots enriched for the current cell type, skip to next
    if (1 %notin% Idents(s)) {
      print(paste0("No ", type))
      next
    }
    
    # Subset for current cell type
    type_s <- subset(s, idents = 1)
    
    # Extract raw counts
    counts <- GetAssayData(type_s, assay = "Spatial", layer = "counts")
    
    # If comparing iAD-Extensive or iAD-Limited, set idents to condition_clearance
    if ((ident.1 %in% c("ext", "lim")) | (ident.2 %in% c("ext", "lim"))) {
      Idents(type_s) <- "condition_clearance"
    } else {
      Idents(type_s) <- "condition"
    }
    
    # Calculate percent expression for all genes based on raw counts  
    LFC <- FoldChange(object = type_s, ident.1 = ident.1, ident.2 = ident.2, fc.name = "avg_log2FC", base = 2, layer = "counts") %>% data.frame()
    
    # Filter for genes expressed in at least 1% of either group
    genes_keep <- row.names(LFC)[LFC$pct.1 >= 0.01 | LFC$pct.2 >= 0.01]
    
    # Remove contamination
    if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
      genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)] 
    }
    
    # Filter count matrix
    type_counts[[type]] <- counts[row.names(counts) %in% genes_keep,]
    print(dim(type_counts[[type]]))
  }
  comparison_counts[[comparison]] <- type_counts
}

# Extract metadata for each cell type, and store in a list (meta data will be the same for all comparisons)
type_meta <- list()
for (type in cell_types) {
  type_meta[[type]] <- s@meta.data[s@meta.data[[type]] == 1,]
}

# Confirm raw count columns and metadata rows match for each comparison and each cell type
for (comp in comparisons) {
  type_counts <- comparison_counts[[comp]]
  for (type in names(type_counts)) {
    print(sum(colnames(type_counts[[type]]) != row.names(type_meta[[type]])))
  }
}

# Sum raw counts by sample ID to generate pseudobulk data
comparison_bulk_counts <- list()
for (comp in comparisons) {
  type_counts <- comparison_counts[[comp]]
  bulk_counts <- list()
  for (type in names(type_counts)) {
    colnames(type_counts[[type]]) <- type_meta[[type]]$sample_id
    bulk <- t(rowsum(t(type_counts[[type]]), group = colnames(type_counts[[type]])))
    bulk_counts[[type]] <- bulk
  }
  comparison_bulk_counts[[comp]] <- bulk_counts
}

# Subset metadata to have one row per sample, and standardize continuous covariates 
for (type in cell_types) {
  meta <- type_meta[[type]][,c('age', 'sex', 'sample_id', 'gw_avg_features', 'gDNA_percent', 'condition', 'condition_clearance')]
  meta <- unique(meta)
  meta$age_centered <- as.numeric((meta$age - mean(meta$age)) / sd(meta$age))
  meta$gw_features_centered <- as.numeric((meta$gw_avg_features - mean(meta$gw_avg_features)) / sd(meta$gw_avg_features))
  meta$gDNA_centered <- as.numeric((meta$gDNA_percent - mean(meta$gDNA_percent)) / sd(meta$gDNA_percent))
  meta$sex <- as.factor(meta$sex)
  row.names(meta) <- meta$sample_id
  type_meta[[type]] <- meta
}

# Confirm pseudobulk columns and metadata rows match for each comparison and each cell type
for (comp in comparisons) {
  bulk_counts <- comparison_bulk_counts[[comp]]
  for (type in names(bulk_counts)) {
    meta <- type_meta[[type]]
    bulk <- bulk_counts[[type]]
    bulk <- bulk[, row.names(meta)]
    
    print(sum(colnames(bulk) != row.names(meta)))
    print(sum(is.na(bulk)))
    
    bulk_counts[[type]] <- bulk
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
  
  # Run DESeq2 for each cell type
  type_results <- list()
  for (type in names(bulk_counts)) {
    meta <- type_meta[[type]]
    
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
        dds <- DESeqDataSetFromMatrix(countData = bulk_counts[[type]], colData = meta, design= ~ sex + age_centered + gw_features_centered +
                                        gDNA_centered + condition)
        results <- DESeq(dds)
      },
      error = function(e) {
        skip_to_next <- TRUE
      }
    )
    if (skip_to_next) {
      print(paste0("Error in ", comparison, ", ", type))
      next
    }
    
    tryCatch(
      {
        results <- results(results, name = paste0(comp_name))
        type_results[[type]] <- results
      },
      error = function(e) {
        print(paste0("No ", comparison, " results for ", type))
      }
    )
  }
  comparison_results[[comparison]] <- type_results
}

# Set DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Define DEGs
for (comparison in comparisons) {
  type_results <- comparison_results[[comparison]]
  for (type in names(type_results)) {
    print(type)
    results <- type_results[[type]]
    results$gene <- row.names(results)
    results$DE <- "Not DE"
    results$DE[results$log2FoldChange > fc_thresh & results$padj < p_thresh] <- "Upregulated"
    results$DE[results$log2FoldChange < -fc_thresh & results$padj < p_thresh] <- "Downregulated"
    results$DE_gene <- NA
    results$DE_gene[results$DE != "Not DE"] <- results$gene[results$DE != "Not DE"]
    type_results[[type]] <- results
  }
  comparison_results[[comparison]] <- type_results
}

# Save results
for (comparison in comparisons) {
  type_results <- comparison_results[[comparison]]
  for (type in names(type_results)) {
    write.csv(type_results[[type]], paste0(output_folder, comparison, "/deseq/results/", type, ".csv"))
  }
}






    
  


  


