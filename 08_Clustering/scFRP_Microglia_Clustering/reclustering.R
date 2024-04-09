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
# Summary: Clean microglia clustering
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

# Define filter operator 
`%notin%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
# Re-integrate microglia data without clusters 1 and 8

# Load initial microglia cluster object
load("/path/to/s_microglia_new")

# Remove clusters 1 and 8 from original microglia clustering, resolution = 0.45 (cluster 1 similar to neurons, cluster 8 similar to astrocytes based on marker genes)
s <- subset(s, SCT_snn_res.0.45 %notin% c(1, 8))

# Merge data for sample B1 across pools (low cell count after filtering)
s@meta.data$new_layer_id <- s@meta.data$unique_id
s@meta.data$new_layer_id[s@meta.data$sample_id == "NMA22.205.B1"] <- "NMA22.205.B1"
DefaultAssay(s) <- "RNA"
print(table(s$new_layer_id))

# Join layers and then split based on merged layer ID (layers are split by pool and sample ID, except for sample B1) 
s[['RNA']] <- JoinLayers(s[['RNA']])
s[['RNA']] <- split(s[['RNA']], f = s$new_layer_id)

# Rerun SCTransform and PCA using the same parameters as original integration
s <- SCTransform(s, assay = "RNA", ncells = 10000, vars.to.regress = c("percent.mt"))
s <- RunPCA(s, npcs = 30, verbose = F)

# Integrate data (k.weight adjusted to be less than ncells in smallest layer)
s <- IntegrateLayers(object = s, method = CCAIntegration, normalization.method = 'SCT', orig.reduction = 'pca', verbose = F, k.weight = 40)

# Re-generate UMAP (use 20 PCs to match prior UMAP dims used)
s <- RunUMAP(s, dims = 1:20, reduction = "integrated.dr", verbose = FALSE, reduction.name = "integratedUMAP", min.dist = 0.1)

#-------------------------------------------------------------------------------
# Initial clustering using a series of resolutions

# Find neighbors (using same number of dims as UMAP)
s <- FindNeighbors(s, reduction = "integrated.dr", dims = 1:20)

# Find clusters using a sequence of resolutions 
for (i in seq(0.05, 1, 0.05)) {
  print(i)
  s <- FindClusters(s, resolution = i, graph.name = 'SCT_snn')
}

# Create variable for clusters using resolution = 0.55 (based on visual analysis of UMAP clusters)
s$microglia_cluster <- paste0("microglia_", s$SCT_snn_res.0.55)
s$microglia_cluster <- factor(s$microglia_cluster, levels = paste0("microglia_", 0:12))

# Save new Seurat object
save(s, file = "/path/to/s_microglia_new_clean")

# Save annotations (for updating subtypes in label transfer object)
write.csv(data.frame(row_name = rownames(s@meta.data), cluster = s$microglia_cluster), 
          paste0(output_dir, "microglia_clusters_new_clean/annotations.csv"))

#-------------------------------------------------------------------------------
# Find marker genes for clusters

# Define output folder, set idents, recorrect SCT counts across layers
marker_dir <- paste0(output_dir, "microglia_clusters_new_clean/markers/")
s <- PrepSCTFindMarkers(s)
Idents(s) <- "microglia_cluster"

# Run FindMarkers
markers <- list()
sig_markers <- list()
for (clust in unique(s@meta.data$microglia_cluster)){
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

