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
# Summary: Initial microglia clustering
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
  library("ggdendroplot")
})

# Define output folder
output_dir <- "/path/to/general/output/folder/"

#-------------------------------------------------------------------------------
# Re-integrate data for cells defined as microglia via label transfer

# Load label transfer Seurat object 
load("/path/to/s_label_transfer_final")

# Subset for microglia subtypes
s <- subset(s, subtype %in% c("Microglia", "Dividing Microglia"))

# Merge data for sample B1 across pools (low cell count after filtering)
s@meta.data$new_layer_id <- s@meta.data$unique_id
s@meta.data$new_layer_id[s@meta.data$sample_id == "NMA22.205.B1"] <- "NMA22.205.B1"
DefaultAssay(s) <- "RNA"

# Join layers and then split based on merged layer ID (layers are split by pool and sample ID, except for sample B1) 
s[['RNA']] <- JoinLayers(s[['RNA']])
s[['RNA']] <- split(s[['RNA']], f = s$new_layer_id)

# Rerun SCTransform and PCA using the same parameters as original integration
s <- SCTransform(s, assay = "RNA", ncells = 10000, vars.to.regress = c("percent.mt"))
s <- RunPCA(s, npcs = 30, verbose = F)

# Integrate data (k.weight adjusted to be less than ncells in smallest layer)
s <- IntegrateLayers(object = s, method = CCAIntegration, normalization.method = 'SCT', orig.reduction = 'pca', verbose = F, k.weight = 50)

# Re-generate integrated UMAP (use 20 PCs to match prior UMAP dims used)
s <- RunUMAP(s, dims = 1:20, reduction = "integrated.dr", verbose = FALSE, reduction.name = "integratedUMAP")

#-------------------------------------------------------------------------------
# Initial clustering using a series of resolutions

# Find neighbors (using same number of dims as UMAP)
s <- FindNeighbors(s, reduction = "integrated.dr", dims = 1:20)

# Find clusters using a sequence of resolutions 
for (i in seq(0.05, 1, 0.05)) {
  print(i)
  s <- FindClusters(s, resolution = i, graph.name = 'SCT_snn')
}

# Save Seurat object
save(s, file = "/path/to/s_microglia_new")

#-------------------------------------------------------------------------------
# Find marker genes for clusters using resolution = 0.45 (based on visual analysis of UMAP clusters)

# Define output folder
marker_dir <- paste0(output_dir, "microglia_clusters_new/markers_0.45/")

# Define temporary cluster variable, set idents, recorrect SCT counts across layers
s$clust_temp <- s$SCT_snn_res.0.45
Idents(s) <- s$clust_temp
s <- PrepSCTFindMarkers(s)

# Run FindMarkers
markers <- list()
sig_markers <- list()
for (clust in unique(s@meta.data$clust_temp)){
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

# Define output file names
sig_out_file <- paste0(marker_dir, "all_sig_markers.csv")
out_file <- paste0(marker_dir, "all_markers.csv")

# Remove MT and HB genes
all_sig_markers_filtered <- all_sig_markers[!grepl("^MT-", all_sig_markers$gene),]
all_sig_markers_filtered <- all_sig_markers_filtered[!grepl("^HB", all_sig_markers_filtered$gene),]

# Save results
write.csv(all_sig_markers_filtered, sig_out_file, row.names = FALSE, quote = FALSE)
write.csv(all_markers, out_file, row.names = FALSE, quote = FALSE)




















