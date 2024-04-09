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
# Summary: Label transfer for vascular and glial subtypes
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
# Prepare data for subtype label transfer

# Load integrated Seurat object and set default assay
load("/path/to/s_integrated_v5")
DefaultAssay(s) <- "RNA"

# Load broad cell type annotations
anno <- read.csv(paste0(output_dir, "label_transfer_new/annotations/broad_annotations.csv"), row.names = 1)

# Load ROSMAP reference object (5000 cells per broad cell type, same reference used for broad annotations)
ref <- readRDS("/path/to/ROSMAP/general/reference/GW.rds")

# Initialize a clean Seurat object with the integrated data and meta data
counts <- GetAssayData(ref, assay = "RNA", layer = "counts")
data <- GetAssayData(ref, assay = "RNA", layer = "data")
new_ref <- CreateSeuratObject(counts = counts, data = data)
new_ref@meta.data <- ref@meta.data

# Save memory 
ref <- NULL
gc()

#-------------------------------------------------------------------------------
# Annotate glial subtypes (query object layers are split by pool and sample ID)

# Subset for glial cell types (with broad prediction score = 1)
temp_s <- subset(s, cells = row.names(anno)[anno$prediction_score == 1 & anno$predicted_type %in% c("Oligodendrocyte", "OPC", "Astrocyte")])
temp_ref <- subset(new_ref, CellTypeBroad %in% c("Oligodendrocyte", "OPC", "Astrocyte"))

# Find 3000 variable features and run PCA on the subset reference
temp_ref <- temp_ref %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA()
print(sum(VariableFeatures(temp_ref) %in% rownames(temp_s)))

# Find transfer anchors (query object layers are split by sample/pool)
anchors <- FindTransferAnchors(reference = temp_ref, query = temp_s, normalization.method = "LogNormalize", 
                               query.assay = "RNA", reference.assay = "RNA", reference.reduction = "pca", dims = 1:30)

# Save anchors
saveRDS(anchors, paste0(output_dir, "label_transfer_new/annotations/glia_anchors.rds"))

# Transfer data for glial subtypes
temp_s <- TransferData(anchorset = anchors, reference = temp_ref, query = temp_s, refdata = temp_ref$ct, k.weight = 20)

# Save cell type predictions and scores
write.csv(data.frame(sample_id = temp_s$sample_id, barcode = str_split_fixed(row.names(temp_s@meta.data), "_", 3)[,3], 
                     predicted_type = temp_s$predicted.id, prediction_score = temp_s$predicted.id.score, row.names = row.names(temp_s@meta.data)), 
          paste0(output_dir, "label_transfer_new/annotations/glia_annotations.csv"))

#-------------------------------------------------------------------------------
# Annotate vascular subtypes (query object layers are split by pool and sample ID)

# Subset for vascular cell types (with broad prediction score = 1)
temp_s <- subset(s, cells = row.names(anno)[anno$prediction_score == 1 & anno$predicted_type %in% c("Endothelial", "Stromal")])
temp_ref <- subset(new_ref, CellTypeBroad %in% c("Endothelial", "Stromal"))

# Find 3000 variable features and run PCA on the subset reference
temp_ref <- temp_ref %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA()
print(sum(VariableFeatures(temp_ref) %in% rownames(temp_s)))

# Find transfer anchors (query object layers are split by sample/pool)
anchors <- FindTransferAnchors(reference = temp_ref, query = temp_s, normalization.method = "LogNormalize", 
                               query.assay = "RNA", reference.assay = "RNA", reference.reduction = "pca", dims = 1:30)

# Save anchors
saveRDS(anchors, paste0(output_dir, "label_transfer_new/annotations/vasculature_anchors.rds"))

# Transfer data for vascular subtypes
temp_s <- TransferData(anchorset = anchors, reference = temp_ref, query = temp_s, refdata = temp_ref$ct, k.weight = 20)

# Save cell type predictions and scores
write.csv(data.frame(sample_id = temp_s$sample_id, barcode = str_split_fixed(row.names(temp_s@meta.data), "_", 3)[,3], 
                     predicted_type = temp_s$predicted.id, prediction_score = temp_s$predicted.id.score, row.names = row.names(temp_s@meta.data)), 
          paste0(output_dir, "label_transfer_new/annotations/vasculature_annotations.csv"))










