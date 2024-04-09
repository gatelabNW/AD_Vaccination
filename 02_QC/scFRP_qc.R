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
# Summary: QC with DoubletFinder and Seurat
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("ggthemes")
  library("ggrepel")
  library("grid")
  library("DoubletFinder")
  library("doMC")
  library("xlsx")
  library("RColorBrewer")
})

# Define function to run DoubletFinder
run_doubletfinder <- function(s) {
  
  # Standard data normalization and scaling
  s <- s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
  
  # Run PCA, generate basic clusters, run TSNE
  # Note that clustering is not used since we are not adjusting for homotypic doublet proportion
  s <- RunPCA(s, features = VariableFeatures(object = s))
  s <- FindNeighbors(s, dims = 1:10)
  s <- FindClusters(s) 
  s <- RunTSNE(s, dims = 1:10)
  
  # pK Identification (not using ground-truth classifications)
  sweep.res.list <- paramSweep(s, PCs = 1:10, sct = FALSE) 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))
  
  # Homotypic Doublet Proportion Estimate (ultimately not used to identify doublets)
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_formation_rate*nrow(s@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # Run DoubletFinder  
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # Rename meta data column for consistency
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.classifications"
  print(table(s[["DF.classifications"]]))
  
  return(s)
}

# Define doublet formation rate
# Average cells per probe barcode is 2813 (manually calculated from CellRanger output htmls)
# Estimated doublet formation rate is 2.25% (scaled from estimated undetectable multiplet rate of 0.4% for 500 cells recovered/probe barcode)
doublet_formation_rate <- 0.0225

#-------------------------------------------------------------------------------
# Initialize Seurat objects and run DoubletFinder on cohort 5 scFRP pools/samples

# Generate raw merged Seurat objects for pool 1 and pool 2
for (pool in c(1, 2)) {
  
  # Define path to SoupX output for current pool
  soupx_dir <- paste0("/path/to/SoupX/output/folder/pool_", pool, "/corrected_counts")

  # Create output folder for current pool
  output_dir <- paste0("/path/to/QC/output/pool_", pool, "/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define paths to SoupX corrected counts for each sample in the current pool
  soupx_sample_dirs <- list.dirs(path = soupx_dir, full.names = TRUE, recursive = FALSE)
  samples <- list.dirs(path = soupx_dir, full.names = FALSE, recursive = FALSE)
  samples <- str_replace_all(samples, "-", ".")
  print(samples)
  
  # Define function to initialize Seurat objects
  load_seurat <- function(dir) {
    counts <- Read10X(dir)
    project <- tail(unlist(strsplit(dir, "/")), 1)
    project <- str_replace_all(project, "-", ".")
    return(CreateSeuratObject(counts = counts, project = project)) 
  }
  
  # Initialize Seurat objects for all samples, using SoupX corrected counts
  seurat_object_list <- sapply(soupx_sample_dirs, load_seurat)
  
  # Run DoubletFinder on all samples
  seurat_object_list <- lapply(seurat_object_list, run_doubletfinder)
  
  # Merge objects and add sample ID to cell names
  s <- merge(seurat_object_list[[1]],
             unlist(seurat_object_list,
                    use.names = FALSE)[2:length(seurat_object_list)],
             add.cell.ids = samples,
             project = "LCMB")
  
  # Save merged Seurat object for current pool 
  save(s, file = paste0(output_dir, "s_raw_merged"))
}

#-------------------------------------------------------------------------------
# Filter cells based on percent mitochondrial expression and remove doublets

# QC filtering for each sample in each pool
for (pool in c(1, 2)) {
  
  # Define output folder for current pool
  output_dir <- paste0("/path/to/QC/output/pool_", pool, "/")
  
  # Load raw merged object for current pool
  load(file = paste0(output_dir, "s_raw_merged"))
  
  # Calculate Percent Mitochondrial expression
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
  
  # Cells per sample pre doublet filtering
  pre_doublet <- table(s$orig.ident)
  
  # Remove doublets
  s <- subset(s, subset = DF.classifications == "Singlet")
  
  # Cells per sample post doublet filtering
  post_doublet <- table(s$orig.ident)
  
  # Table of doublets removed per sample
  doublets_removed <- pre_doublet - post_doublet
  
  # Percent doublets removed
  pct_doublets_removed <- (doublets_removed / pre_doublet)*100
  
  # Save doublet QC stats
  write.csv(doublets_removed, paste0(output_dir, "doublets_removed.csv"))
  write.csv(pct_doublets_removed, paste0(output_dir, "pct_doublets_removed.csv"))
  
  # Cells per sample pre MT filtering
  pre_mt <- table(s$orig.ident)
  
  # Filter based on MT % 
  s <- subset(s, subset = percent.mt < 20)
  
  # Cells per sample post MT filtering
  post_mt <- table(s$orig.ident)
  
  # Table of MT spots removed per sample
  mt_removed <- pre_mt - post_mt
  
  # Percent MT spots removed
  pct_mt_removed <- (mt_removed / pre_mt)*100
  
  # Save MT QC stats
  write.csv(mt_removed, paste0(output_dir, "mt_spots_removed.csv"))
  write.csv(pct_mt_removed, paste0(output_dir, "pct_mt_spots_removed.csv"))
  
  # Save Seurat object
  save(s, file = paste0(output_dir, "s_post_qc_merged"))
}



