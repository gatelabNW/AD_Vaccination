# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                         AN1792 Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 04-04-2024
# Written by: Anne Forsyth
# Summary: Identify markers for cortical amyloid clusters
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
  library("randomcoloR")
  library("ggpubr")
  library("gridExtra")
  library("ggplot2")
  library("ggnewscale")
})

# Define output folder
output_dir <- "/path/to/findmarkers/results/"

# Load Seurat object with WNN
s <- readRDS("/path/to/plaque_wnn.rds")

#-------------------------------------------------------------------------------
# Run FindMarkers for WNN clusters, resolution = 0.5

# Name of meta data column with cluster annotations
cluster_name <- "wsnn_res.0.5"

# Set Idents 
Idents(s) <- cluster_name

# Re-correct SCT counts across samples 
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}
s <- PrepSCTFindMarkers(s)

# Initialize lists for cluster results
markers <- list()
sig_markers <- list()

# Run FindMarkers on each cluster 
for (clust in unique(s@meta.data[,cluster_name])){
  print(clust)
  
  results <- FindMarkers(s, ident.1 = clust, min.pct = 0.25, recorrect_umi = FALSE, only.pos = TRUE)
  results$BH <- p.adjust(results$p_val, method = "BH")
  results$cluster <- clust
  results$gene <- rownames(results)
  results$PFC <- -log10(results$BH)*abs(results$avg_log2FC)
  sig_results <- dplyr::filter(results, BH < 0.05)
  
  markers[[clust]] <- results
  sig_markers[[clust]] <- sig_results
}

# Compile results
all_sig_markers <- data.table::rbindlist(sig_markers)|>as.data.frame()
all_markers <- data.table::rbindlist(markers)|>as.data.frame()

# Define output file paths
sig_out_file <- paste0(output_dir, "all_sig_markers.csv")
out_file <- paste0(output_dir, "all_markers.csv")

# Remove contamination
all_sig_markers_filtered <- all_sig_markers[!grepl("^MT-", all_sig_markers$gene),]
all_sig_markers_filtered <- all_sig_markers_filtered[!grepl("^HB", all_sig_markers_filtered$gene),]

# Save results
write.csv(all_sig_markers_filtered, sig_out_file, row.names = FALSE, quote = FALSE)
write.csv(all_markers, out_file, row.names = FALSE, quote = FALSE)

#-------------------------------------------------------------------------------
# Find positive and negative markers for all clusters (used for Fig 4L)

# Run FindMarkers for each cluster, not filtering results by adjusted p-value
markers <- list()
for (clust in unique(s@meta.data$cluster)){
  print(clust)
  
  results <- FindMarkers(s, ident.1 = clust, min.pct = 0.25, recorrect_umi = FALSE, only.pos = FALSE)
  
  results$BH <- p.adjust(results$p_val, method = "BH")
  results$cluster <- clust
  results$gene <- rownames(results)
  results$PFC <- -log10(results$BH)*abs(results$avg_log2FC)
  
  markers[[clust]] <- results
}

# Compile results
all_markers <- data.table::rbindlist(markers)|>as.data.frame()

# Remove contamination
filtered <- all_markers[!grepl("^MT-", all_markers$gene),]
filtered <- filtered[!grepl("^HB", filtered$gene),]

# Save results
write.csv(filtered, paste0(output_dir, "pos_neg_markers_unfiltered/all_markers.csv"), row.names = FALSE, quote = FALSE)

