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
# Summary: RNA integration with Seurat
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
})

# Define input/output folders
prev_out_dir <- "/path/to/qc/output/folder"
objects_out_dir <- "/path/to/data/output/folder" 

# Load RNA and protein post-QC objects
all_seurat <- glue("{prev_out_dir}/all_samples_02_rna.rds")|>readRDS()
pro_filtered <- readRDS(paste0(prev_out_dir, "/all_samples_02_pro_filtered.rds"))

# Define filter operator
`%notin%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
# Integrate cohort 5/7 RNA data (layers are split by sample ID)

# Remove spots not located on cortical or hippocampal tissue
spots_remove_1 <- all_seurat@meta.data$sample_barcode[all_seurat@meta.data$sample_id %in% c("NMA22.A9", "NMA22.B9") & 
                                                            all_seurat@meta.data$Manual_Layer %in% c("exclude", "non-hippocampus")]
spots_remove_2 <- all_seurat@meta.data$sample_barcode[all_seurat@meta.data$sample_id %notin% c("NMA22.A9", "NMA22.B9") & 
                                                            all_seurat@meta.data$Manual_Layer == "exclude"]
spots_remove <- c(spots_remove_1, spots_remove_2)

all_seurat <- subset(all_seurat, sample_barcode %notin% spots_remove)
pro_filtered <- subset(pro_filtered, sample_barcode %notin% spots_remove)

# Create temporary objects with merged sample data to identify spots with zero protein or gene expression
temp_rna <- JoinLayers(all_seurat)
temp_pro_filtered <- JoinLayers(pro_filtered)

# There are no spots with zero RNA expression
spots.remove <- colnames(temp_rna[, colSums(temp_rna@assays$Spatial@layers$counts)==0])

# Remove spots with zero protein expression 
spots.remove <- colnames(temp_pro_filtered[, colSums(temp_pro_filtered@assays$Protein@layers$counts)==0])
all_seurat <- subset(all_seurat, sample_barcode %notin% spots.remove)
pro_filtered <- subset(pro_filtered, sample_barcode %notin% spots.remove)

# Save list of spots removed based on zero protein expression
write.csv(data.frame(spot = spots.remove), paste0(objects_out_dir, "/zero_protein_expression.csv"), row.names = FALSE)

# Save memory
temp_rna <- NULL
temo_pro_filtered <- NULL
gc()

# Save protein object
saveRDS(pro_filtered, paste0(objects_out_dir, "/all_samples_03_pro_filtered.rds"))

# Run SCTransform
all_seurat <- SCTransform(all_seurat, assay = 'Spatial')

# Run PCA
all_seurat <- RunPCA(all_seurat, npcs = 30, verbose = F)

# Integrate data 
all_seurat <- IntegrateLayers(object = all_seurat, method = CCAIntegration, normalization.method = 'SCT', orig.reduction = 'pca', verbose = F)

# Generate integrated UMAP
all_seurat <- RunUMAP(all_seurat, reduction = 'integrated.dr', dims = 1:30, verbose = FALSE, reduction.name = 'integratedUMAP')

# Save object
saveRDS(all_seurat, paste0(objects_out_dir, "/all_samples_03_rna.rds"))





