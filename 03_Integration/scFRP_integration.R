# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AN1792 Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 04-03-2024
# Written by: Anne Forsyth
# Summary: Integration with Seurat
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
})

# Define output folder
output_dir <- "/path/to/output/folder/"

# Load pool 1 and pool 2 Seurat objects
load("/path/to/pool1/qc/merged/object")
s1 <- s
load("/path/to/pool2/qc/merged/object")
s2 <- s

# Save memory
s <- NULL
gc()

#-------------------------------------------------------------------------------
# Integrate cohort 5 scFRP data (layers are split by pool and sample ID)

# Rename cells so merged cell names will be unique 
s1 <- RenameCells(s1, add.cell.id = "pool1")
s2 <- RenameCells(s2, add.cell.id = "pool2")

# Merge pool 1 and pool 2 objects
s <- merge(s1, s2, merge.data = TRUE)

# Save memory
s1 <- NULL
s2 <- NULL
gc()

# Add meta data columns for pool, sample ID, unique ID 
s@meta.data$pool <- str_split_fixed(row.names(s@meta.data), "_", 2)[,1]
s@meta.data$sample_id <- s@meta.data$orig.ident
s@meta.data$unique_id <- paste0(s@meta.data$pool, "_", s@meta.data$sample_id)

# Run SCTransform 
s <- SCTransform(s, assay = "RNA", ncells = 10000, vars.to.regress = c("percent.mt"))

# Run PCA
s <- RunPCA(s, npcs = 30, verbose = F)

# Integrate data 
s <- IntegrateLayers(object = s, method = CCAIntegration, normalization.method = 'SCT', orig.reduction = 'pca', verbose = F)

# Elbow plot to determine number of PCs to use 
ElbowPlot(s, ndims = 30)

# Generate integrated UMAP
s <- RunUMAP(s, dims = 1:20, reduction = "integrated.dr", verbose = FALSE, reduction.name = "integratedUMAP")

# Save integrated object
save(s, file = paste0(output_dir, "s_integrated_v5"))


