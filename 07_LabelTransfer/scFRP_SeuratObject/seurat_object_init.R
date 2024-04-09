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
# Summary: Initialize Seurat object with broad label transfer annotations
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
})

# Define output folder
output_dir <- "/path/to/general/output/folder/"

# Load integrated Seurat object
load("/path/to/s_integrated_v5")

# Remove cluster 14 (similar to red blood cells, not in reference) and set default assay 
s <- subset(s, SCT_snn_res.0.5 != 14)
DefaultAssay(s) <- "RNA"

# Load ROSMAP reference object used for broad label transfer
ref <- readRDS("/path/to/ROSMAP/general/reference/GW.rds") 

#-------------------------------------------------------------------------------
# Prepare the reference the same as broad label transfer and project query cells onto broad reference UMAP

# Initialize a clean Seurat object with the integrated data and meta data
counts <- GetAssayData(ref, assay = "RNA", layer = "counts")
data <- GetAssayData(ref, assay = "RNA", layer = "data")
new_ref <- CreateSeuratObject(counts = counts, data = data)
new_ref@meta.data <- ref@meta.data

# Save memory 
ref <- NULL
gc()

# Find 3000 variable features to increase overlap with query dataset
new_ref <- FindVariableFeatures(new_ref, nfeatures = 3000)
print(sum(VariableFeatures(new_ref) %in% rownames(s)))

# Run PCA using 3000 variable features 
new_ref <- new_ref %>% ScaleData() %>% RunPCA()

# Load broad cell type anchors and annotations
anchors <- readRDS(paste0(output_dir, "label_transfer_new/annotations/broad_anchors.rds"))
anno <- read.csv(paste0(output_dir, "label_transfer_new/annotations/broad_annotations.csv"), row.names = 1)
sum(rownames(anno) != rownames(s@meta.data))

# Add annotations to seurat object
s$broad_type <- anno$predicted_type
s$broad_score <- anno$prediction_score

# Generate reference UMAP and project query cells 
new_ref <- RunUMAP(new_ref, dims = 1:30, reduction = 'pca', return.model = TRUE)
s <- IntegrateEmbeddings(anchorset = anchors, reference = new_ref, query = s, new.reduction.name = "ref.pca")
s <- ProjectUMAP(query = s, query.reduction = "ref.pca", reference = new_ref, reference.reduction = "pca", reduction.model = "umap")

#-------------------------------------------------------------------------------
# Label "unknown" cells, add subtype annotations, generate UMAP

# Label broad "Unknown" cells
s$broad_type[s$broad_score < 1] <- "Unknown"

# Regenerate integrated UMAP for cells used for label transfer
s <- RunUMAP(s, dims = 1:20, reduction = "integrated.dr", verbose = FALSE, reduction.name = "integratedUMAP")

# Load immune annotations
immune_anno <- read.csv(paste0(output_dir, "label_transfer_new/annotations/immune_ct_annotations.csv"), row.names = 1)
tnk <- read.csv(paste0(output_dir, "label_transfer_new/annotations/tnk_immune_annotations.csv"), row.names = 1)
tnk <- tnk[tnk$predicted_type == "T/NK",]
immune_anno <- rbind(immune_anno, tnk)

# Load glial and vascular annotations
annotations <- data.frame()
for (type in c("glia", "vasculature")) {
  anno <- read.csv(paste0(output_dir, "label_transfer_new/annotations/", type, "_annotations.csv"), row.names = 1)
  annotations <- rbind(annotations, anno)
}

# Combine all cell type annotations
annotations <- rbind(annotations, immune_anno)
length(unique(row.names(annotations))) == nrow(annotations)

# Add "unknown" labels for unknown broad cell types, add neuronal subtypes
cells_add <- row.names(s@meta.data)[s@meta.data$broad_type == "Unknown" | s@meta.data$broad_type %in% c("EN", "Interneuron")]
add_meta <- s@meta.data[cells_add,]
add_meta$broad_score[add_meta$broad_type == "Unknown"] <- NA # Subtype score is NA for unknown broad cell types

# Combine subtype labels to add
anno_add <- rbind(data.frame(type = annotations$predicted_type, score = annotations$prediction_score, row.names = row.names(annotations)), 
                  data.frame(type = add_meta$broad_type, score = add_meta$broad_score, row.names = rownames(add_meta)))

# Confirm NA are expected
sum(is.na(anno_add)) == sum(s$broad_type == "Unknown")

# Add to meta data
s$subtype <- anno_add$type
s$subtype_score <- anno_add$score

# Annotate subtypes as "Unknown" if subtype prediction score <= 0.5 (interneuron and EN annotations will not change since using broad annotations which have been filtered for score = 1)
s$subtype[s$subtype_score <= 0.5] <- "Unknown"

# Save object
save(s, file = "/path/to/s_label_transfer_final")

