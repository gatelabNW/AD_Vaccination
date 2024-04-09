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
# Summary: Label transfer for immune subtypes
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

# Load integrated Seurat object and set default assay
load("/path/to/s_integrated_v5")
DefaultAssay(s) <- "RNA"

# Load broad cell type annotations
anno <- read.csv(paste0(output_dir, "label_transfer_new/annotations/broad_annotations.csv"), row.names = 1)

#-------------------------------------------------------------------------------
# Prepare reference for immune label transfer with T/NK cells

# Load ROSMAP immune reference with T/NK cells 
ref <- readRDS("/path/to/ROSMAP/tnk/immune_cell_annotation.rds") 

# Initialize a clean Seurat object with the integrated data and meta data ("data" slot contains integrated data, "counts" slot is empty)
data <- GetAssayData(ref, assay = "RNA", layer = "data")
new_ref <- CreateSeuratObject(counts = data, data = data)
new_ref@meta.data <- ref@meta.data

# Find 3000 variable features in the reference to increase overlap with query dataset, run PCA
new_ref <- new_ref %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA()
print(sum(VariableFeatures(new_ref) %in% rownames(s)))

# Save memory 
ref <- NULL
gc()

#-------------------------------------------------------------------------------
# Annotate T/NK cells

# Subset for immune cells (with broad prediction score = 1)
temp_s <- subset(s, cells = row.names(anno)[anno$prediction_score == 1 & anno$predicted_type == "Immune"])

# Find transfer anchors (query object layers are split by sample/pool)
anchors <- FindTransferAnchors(reference = new_ref, query = temp_s, normalization.method = "LogNormalize", 
                               query.assay = "RNA", reference.assay = "RNA", reference.reduction = "pca", dims = 1:30)

# Save anchors
saveRDS(anchors, paste0(output_dir, "label_transfer_new/annotations/tnk_anchors.rds"))

# Transfer data for immune subtypes
temp_s <- TransferData(anchorset = anchors, reference = new_ref, query = temp_s, refdata = new_ref$ct, k.weight = 20)

# Save cell type predictions and scores
write.csv(data.frame(sample_id = temp_s$sample_id, barcode = str_split_fixed(row.names(temp_s@meta.data), "_", 3)[,3], 
                     predicted_type = temp_s$predicted.id, prediction_score = temp_s$predicted.id.score, row.names = row.names(temp_s@meta.data)), 
          paste0(output_dir, "label_transfer_new/annotations/tnk_immune_annotations.csv"))

#-------------------------------------------------------------------------------
# Prepare reference for remaining immune cells 

# Load ROSMAP subset immune reference
ref <- readRDS("/path/to/ROSMAP/human_Mglia_immune.subset.artefacts.excluded.rds")

# Initialize a clean Seurat object with the integrated data and meta data
counts <- GetAssayData(ref, assay = "RNA", layer = "counts")
data <- GetAssayData(ref, assay = "RNA", layer = "data")
new_ref <- CreateSeuratObject(counts = counts, data = data)
new_ref@meta.data <- ref@meta.data

# Find 3000 variable features in the reference to increase overlap with query dataset, run PCA
new_ref <- new_ref %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA()
print(sum(VariableFeatures(new_ref) %in% rownames(s)))

# Save memory 
ref <- NULL
gc()

#-------------------------------------------------------------------------------
# Label transfer for remaining immune cells (query object layers are split by pool and sample ID)

# Remove cells previously annotated as T/NK
temp_s <- subset(temp_s, predicted.id != "T/NK")

# Find transfer anchors 
anchors <- FindTransferAnchors(reference = new_ref, query = temp_s, normalization.method = "LogNormalize", 
                               query.assay = "RNA", reference.assay = "RNA", reference.reduction = "pca", dims = 1:30)

# Save anchors 
saveRDS(anchors, paste0(output_dir, "label_transfer_new/annotations/immune_anchors.rds"))

# Transfer data for immune subtypes
temp_s <- TransferData(anchorset = anchors, reference = new_ref, query = temp_s, refdata = new_ref$ct, k.weight = 20)

# Save cell type predictions and scores
write.csv(data.frame(sample_id = temp_s$sample_id, barcode = str_split_fixed(row.names(temp_s@meta.data), "_", 3)[,3], 
                     predicted_type = temp_s$predicted.id, prediction_score = temp_s$predicted.id.score, row.names = row.names(temp_s@meta.data)), 
          paste0(output_dir, "label_transfer_new/annotations/immune_ct_annotations.csv"))









