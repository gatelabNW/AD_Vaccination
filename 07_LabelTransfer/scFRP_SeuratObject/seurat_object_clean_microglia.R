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
# Summary: Refine subtype annotations in label transfer Seurat object
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
})

# Define general output folder
output_dir <- "/path/to/general/output/folder/"

# Load label transfer Seurat object 
load("/path/to/s_label_transfer_final")

# Load data containing cells used in final microglia clusters  
anno <- read.csv(paste0(output_dir, "microglia_clusters_new_clean/annotations.csv"), row.names = 1)

# Define filter operator
`%notin%` <- Negate(`%in%`)

# Create clean subtype variable, labeling microglia cells unknown if not included in final microglia clusters
s@meta.data$subtype_clean <- s@meta.data$subtype
s@meta.data$subtype_clean[s$subtype_clean %in% c("Microglia", "Dividing Microglia") & rownames(s@meta.data) %notin% anno$row_name] <- "Unknown"

# Save object
save(s, file = "/path/to/s_label_transfer_final")






