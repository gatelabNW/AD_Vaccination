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
# Summary: Clustering of cortical amyloid-rich spots
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages ({
  library('Seurat')
  library('glue')
  library('dplyr')
  library('stringr')
  library("Matrix")
  library("data.table")
  library('sctransform')
  library('ggplot2')
  library("randomcoloR")
  library("ggpubr")
  library("tidyverse")
  library("introdataviz")
  library("scales")
})

# Load Seurat object with integrated protein data
pro_orig <- readRDS("/path/to/all_samples_03_pro_filtered.rds")

# Load Seurat object with integrated RNA data
rna_orig <- readRDS("/path/to/all_samples_03_rna_updated.rds")

# Define cortical amyloid-rich spots
rna_orig$plaque_rich <- "not_rich"
rna_orig$plaque_rich[rna_orig$plaque_fluo > 153 & rna_orig$vessel_fluo == 0] <- "rich"

#-------------------------------------------------------------------------------
# Cluster cortical amyloid spots

# Subset RNA and protein objects for cortical amyloid-rich spots in gray matter
rna <- subset(rna_orig, plaque_rich == "rich" & gray_matter != "not_gray")
pro <- subset(pro_orig, cells = row.names(rna@meta.data))

# Add protein integrated DR to RNA object
rna@reductions$prointegrated.dr <- pro@reductions$integrated.dr

# Find neighbors based on integrated DR, using all integrated features 
rna <- FindMultiModalNeighbors(rna, reduction.list = list("integrated.dr", "prointegrated.dr"), dims.list = list(1:ncol(rna@reductions$integrated.dr), 1:ncol(rna@reductions$prointegrated.dr)))

# Generate WNN UMAP 
rna <- RunUMAP(rna, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

# Generate clusters for sequence of resolutions
for (i in seq(0.05, 1, 0.05)) {
  print(i)
  rna <- FindClusters(rna, resolution = i, graph.name = 'wsnn', algorithm = 3, verbose = FALSE)
}

# Save object with WNN
saveRDS(rna, "/path/to/plaque_wnn.rds")

