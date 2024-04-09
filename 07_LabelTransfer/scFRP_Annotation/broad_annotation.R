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
# Summary: Label transfer for broad cell types
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
})

# Define output folder
output_dir <- "/path/to/general/output/folder/"

#-------------------------------------------------------------------------------
# Broad cell type label transfer (query object layers are split by pool and sample ID)

# Load integrated Seurat object
load("/path/to/s_integrated_v5")

# Remove cluster 14 (similar to red blood cells, not in reference) and set default assay 
s <- subset(s, SCT_snn_res.0.5 != 14)
DefaultAssay(s) <- "RNA"

# Load ROSMAP reference object (5000 cells per broad cell type)
ref <- readRDS("/path/to/ROSMAP/general/reference/GW.rds") 

# Initialize a new reference object with the integrated data and meta data
counts <- GetAssayData(ref, assay = "RNA", layer = "counts")
data <- GetAssayData(ref, assay = "RNA", layer = "data")
new_ref <- CreateSeuratObject(counts = counts, data = data)
new_ref@meta.data <- ref@meta.data

# Save memory 
ref <- NULL
gc()

# Find 3000 variable features in the reference to increase overlap with query dataset
new_ref <- FindVariableFeatures(new_ref, nfeatures = 3000)
print(sum(VariableFeatures(new_ref) %in% rownames(s)))

# Run PCA using 3000 variable features 
new_ref <- new_ref %>% ScaleData() %>% RunPCA()

# Find transfer anchors 
anchors <- FindTransferAnchors(reference = new_ref, query = s, normalization.method = "LogNormalize", 
                               query.assay = "RNA", reference.assay = "RNA", reference.reduction = "pca", dims = 1:30)

# Save anchors
saveRDS(anchors, paste0(output_dir, "label_transfer_new/annotations/broad_anchors.rds"))

# Transfer data for broad ROSMAP cell types
s <- TransferData(anchorset = anchors, reference = new_ref, query = s, refdata = new_ref$CellTypeBroad, k.weight = 20)

# Save cell type predictions and scores
write.csv(data.frame(sample_id = s$sample_id, barcode = str_split_fixed(row.names(s@meta.data), "_", 3)[,3], 
                     predicted_type = s$predicted.id, prediction_score = s$predicted.id.score, row.names = row.names(s@meta.data)), 
          paste0(output_dir, "label_transfer_new/annotations/broad_annotations.csv"))





