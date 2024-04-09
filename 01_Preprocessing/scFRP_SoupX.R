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
# Summary: Preprocess scFRP RNA data with SoupX 
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("Matrix")
  library("SoupX")
  library("DropletUtils")
  library("purrr")
  library("Hmisc")
})

# Define function to run SoupX
run_soupx <- function(dir) {
  sample <- unlist(strsplit(dir, "/")) %>%
    tail(1)
  
  sample_dir <- paste0(output_dir, "/corrected_counts/", sample, "/")
  dir.create(sample_dir, showWarnings = FALSE)
  
  print(paste0("Processing sample ", sample))
  
  # Initialize Seurat object for current sample
  seurat <- Read10X(paste0(dir, "/count/sample_filtered_feature_bc_matrix")) %>%
    CreateSeuratObject()
  seurat
  
  # Run SCTransform and PCA
  seurat <- SCTransform(seurat, verbose = TRUE) %>%
    RunPCA()
  
  # Basic clustering
  seurat <- RunUMAP(seurat, dims = 1:10) %>%
    FindNeighbors(dims = 1:10) %>%
    FindClusters() 
  
  # Load 10X data for SoupX 
  toc = Seurat::Read10X(paste0(dir,"/count/sample_filtered_feature_bc_matrix"))
  tod = Seurat::Read10X(paste0(dir,"/count/sample_raw_feature_bc_matrix"))
  
  # Filter table of droplets (tod) so that it has same genes as table of counts (toc)
  tod <- tod[which(rownames(tod) %in% rownames(toc)),]
  
  # Create SoupX object 
  soupx = SoupChannel(tod, toc)
  
  # Add Seurat cluster info to SoupX object
  soupx <- setClusters(soupx, setNames(seurat$seurat_clusters, colnames(seurat)))
  
  # Add Seurat UMAP to SoupX object
  soupx <- setDR(soupx, seurat@reductions$umap[[]], reductName='umap')
  
  # Automatically estimate contamination fraction
  pdf(file = paste0(plot_dir, sample, "_rho_distribution.pdf"), width = 4, height = 4)
  soupx <- autoEstCont(soupx)
  dev.off()
  
  # Print contamination fraction for current sample
  print(paste0("Fraction: ", sample, " ", soupx$metaData$rho[1])) 
  
  # Adjust counts
  adj_counts <- adjustCounts(soupx)
  
  # Save corrected counts
  DropletUtils:::write10xCounts(paste0(output_dir, "/corrected_counts/", sample), adj_counts, overwrite = TRUE)
  
  # Return contamination fraction for current sample
  return(soupx$metaData$rho[1])
}

# Run SoupX for each sample in each pool
for (pool in c(1, 2)) {
  
  # Define path to general CellRanger output for current pool
  cellranger_dir <- paste0("/path/to/cellranger/pool", pool, "/outs/per_sample_outs")
  
  # Create folders for SoupX output and QC plots for current pool
  output_dir <- paste0("/path/to/SoupX/output/folder/pool_", pool)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  plot_dir <- paste0(output_dir, "/qc_plots/")
  dir.create(plot_dir, showWarnings = FALSE)
  
  # Define paths to sample-level CellRanger output for current pool
  sample_dirs <- list.dirs(cellranger_dir, recursive = FALSE)
  print(sample_dirs)
  
  # Run SoupX on all samples in current pool
  rho <- lapply(sample_dirs, run_soupx)
  rho
  
  # Export contamination fraction for each sample (Ext Fig 3B)
  write.csv(rho, file = paste0(output_dir, "/contamination_fraction.csv"))
}

