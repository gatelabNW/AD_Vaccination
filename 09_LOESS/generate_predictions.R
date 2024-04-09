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
# Summary: Generate LOESS predictions for gene expression over amyloid density in gray matter
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggrepel")
  library("ggthemes")
  library("ggpubr")
  library("Seurat")
  library("factoextra")
  library("scales")
  library("patchwork")
  library("clusterProfiler")
  library("pheatmap")
  library("colorspace")
  library("stringi")
  library("randomcoloR")
})

# Define output folder
output_folder <- "/path/to/loess/output/folder/"

# Load microglia gene sets
eggen <- read.csv("/path/to/galatro/et/al/eggen-microglia.csv")
eggen_genes <- eggen$Gene
ham <- read.csv("/path/to/sun/et/al/ham_genes.csv")
ham_genes <- ham$gene
microglia_genes <- unique(c(eggen_genes, ham_genes))

# Load integrated Seurat object 
s <- readRDS("/path/to/all_samples_03.rds")

# Define filter operator 
`%notin%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
# Prepare nAD, iAD-Limited and iAD-Extensive data for LOESS

# Remove amyloid "exclude" spots and subset for gray matter
s <- subset(s, amyloid_filter == "include" & manual_annotation %notin% c("white", "meninges"))

# Recorrect SCT counts across samples (layers are split by sample ID)
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}
s <- PrepSCTFindMarkers(object = s)

# Extract SCT data and meta data
all_data <- GetAssayData(s, assay = "SCT", layer = "data")
all_data <- all_data[row.names(all_data) %in% microglia_genes,]
meta <- s@meta.data

# Save memory 
s <- NULL
gc()

# Convert expression matrix to data frame
data <- data.frame(all_data)

# Save memory 
all_data <- NULL
gc()

# Format data for LOESS
data <- t(data)
data <- data.frame(data) 
data$row_name <- stri_replace_last_fixed(row.names(data), ".", "-")
row.names(data) <- data$row_name
sum(row.names(data) != row.names(meta))
data$amyloid <- meta$amyloid_fluo
data$condition <- meta$condition_clearance

# Define function to generate LOESS predictions 
run_loess <- function(gene, data) {
  
  # Extract gene expression and amyloid density from input data
  data <- data[, c(gene, "amyloid")]
  
  # Generate LOESS model 
  lo <- loess(get(gene) ~ amyloid, data)
  
  # Generate predicted values 
  lo_predict <- predict(lo, data.frame(amyloid = seq(min(data$amyloid), max(data$amyloid), 1))) %>% as.data.frame()
  colnames(lo_predict) <- gene
  
  return(lo_predict)
}

#-------------------------------------------------------------------------------
# Generate LOESS predictions for microglia-associated genes

# Ensure names of genes to test exactly match column names 
microglia_genes <- colnames(data)[colnames(data) %notin% c("row_name", "amyloid", "condition")]

# Generate predictions for iAD-Limited
cur_data <- data[data$condition == "lim", c(microglia_genes, "amyloid")]
predictions <- lapply(microglia_genes, run_loess, data = cur_data)
merged <- as.data.frame(predictions)
merged$amyloid <- seq(min(cur_data$amyloid), max(cur_data$amyloid), 1)
merged <- dplyr::relocate(merged, amyloid) 
write.csv(merged, paste0(output_folder, "lim_predictions.csv"), row.names = FALSE)

# Generate predictions for iAD-Extensive
cur_data <- data[data$condition == "ext", c(microglia_genes, "amyloid")]
predictions <- lapply(microglia_genes, run_loess, data = cur_data)
merged <- as.data.frame(predictions)
merged$amyloid <- seq(min(cur_data$amyloid), max(cur_data$amyloid), 1)
merged <- dplyr::relocate(merged, amyloid) 
write.csv(merged, paste0(output_folder, "ext_predictions.csv"), row.names = FALSE)

# Generate predictions for nAD 
cur_data <- data[data$condition == "nAD", c(microglia_genes, "amyloid")]
predictions <- lapply(microglia_genes, run_loess, data = cur_data)
merged <- as.data.frame(predictions)
merged$amyloid <- seq(min(cur_data$amyloid), max(cur_data$amyloid), 1)
merged <- dplyr::relocate(merged, amyloid) 
write.csv(merged, paste0(output_folder, "nAD_predictions.csv"), row.names = FALSE)

