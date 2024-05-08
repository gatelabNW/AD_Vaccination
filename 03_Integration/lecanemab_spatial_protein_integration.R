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
# Summary: Protein integration with Seurat
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
})

# Define output folder
objects_out_dir <- "/path/to/data/output/folder" 

# Load protein Seurat object with isotype control-normalized data (spots with zero expression have already been removed prior to RNA integration)
pro <- readRDS(paste0(objects_out_dir, "/all_samples_03_pro_filtered.rds"))

#-------------------------------------------------------------------------------
# Integrate cohort 5/7 protein data (layers are split by sample ID)

# CLR-normalize data across cells
pro <- NormalizeData(pro, normalization.method = "CLR", margin = 2)

# Set variable features to all features and scale data
VariableFeatures(pro) <- rownames(pro) 
pro <- ScaleData(pro)

# Generate pseudo-PCs using scaled, CLR-normalized, isotype-normalized data
pseudo <- t(GetAssayData(pro, layer = 'scale.data', assay = 'Protein')) 
colnames(pseudo) <- paste('PC', 1:ncol(pseudo), sep = '_')
pseudo_pca <- CreateDimReducObject(embeddings = as.matrix(pseudo), key = 'pseudoPCA', assay = 'Protein')
pro@reductions$pseudoPCA <- pseudo_pca

# Integrate data using pseudo-PCs and all antibodies
pro <- IntegrateLayers(object = pro, method = RPCAIntegration, assay = "Protein", orig.reduction = "pseudoPCA", verbose = FALSE,
                       features = rownames(pro))

# Generate integrated UMAP
pro <- RunUMAP(pro, reduction = 'integrated.dr', dims = 1:ncol(pro@reductions$integrated.dr), verbose = FALSE, reduction.name = 'integratedUMAP')

# Save integrated object
saveRDS(pro, paste0(objects_out_dir, "/all_samples_03_pro_filtered.rds"))


