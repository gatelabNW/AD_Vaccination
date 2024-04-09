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
# Summary: Initial Seurat clustering to identify cells to exclude from label transfer
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
})

# Define output folder for FindMarkers results
marker_dir <- "/path/to/findmarkers/output/"

# Load integrated Seurat object
load("/path/to/s_integrated_v5")

#-------------------------------------------------------------------------------
# Seurat clustering using a series of resolutions

# Find neighbors (using same dims as integrated UMAP)
s <- FindNeighbors(s, reduction = "integrated.dr", dims = 1:20)

# Find clusters using a sequence of resolutions, print plots to compare resolutions
for (i in seq(0.1, 1, 0.05)) {
  print(i)
  s <- FindClusters(s, resolution = i, graph.name = 'SCT_snn')
  colors <- distinctColorPalette(length(unique(s@meta.data$seurat_clusters)))
  print(DimPlot(s, reduction = "integratedUMAP", group.by = "seurat_clusters", label = TRUE) + scale_color_manual(values = colors) + 
          ggtitle(paste0("Resolution: ", i)))
}

# Save object
save(s, file = "/path/to/s_integrated_v5")  

#-------------------------------------------------------------------------------
# Identify marker genes for clusters using resolution = 0.5

# Format cluster names for FindMarkers
s@meta.data$cluster <- paste0("clust_", s@meta.data$SCT_snn_res.0.5) 

# Set default assay and idents
DefaultAssay(s) <- "SCT"
Idents(s) <- "cluster"

# Recorrect SCT counts across layers
s <- PrepSCTFindMarkers(object = s)

# Find markers for each cluster 
markers <- list()
sig_markers <- list()
for (cluster in unique(s@meta.data$cluster)){
  print(cluster)
  results <- FindMarkers(s, ident.1 = cluster, min.pct = 0.25, recorrect_umi = FALSE, only.pos = TRUE)
  results[["BH"]] <- p.adjust(results$p_val, method = "BH")
  results[["cluster"]] <- cluster
  results[["gene"]] <- rownames(results)
  sig_results <- dplyr::filter(results, BH < 0.05)
  markers[[cluster]] <- results
  sig_markers[[cluster]] <- sig_results
}

# Compile results
all_sig_markers <- data.table::rbindlist(sig_markers)|>as.data.frame()
all_markers <- data.table::rbindlist(markers)|>as.data.frame()

# Define output file paths
sig_out_file <- paste0(marker_dir, "all_sig_markers.csv")
out_file <- paste0(marker_dir, "all_markers.csv")

# Remove contamination
genes_keep <- all_sig_markers$gene
if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
  genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)] 
}
all_sig_markers_filtered <- all_sig_markers[all_sig_markers$gene %in% genes_keep,]

# Save results (all_sig_markers_filtered used to identify top markers)
write.csv(all_sig_markers_filtered, sig_out_file, row.names = FALSE, quote = FALSE)
write.csv(all_markers, out_file, row.names = FALSE, quote = FALSE)








