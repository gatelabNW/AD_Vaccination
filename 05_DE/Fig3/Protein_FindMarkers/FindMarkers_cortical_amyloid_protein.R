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
# Summary: Cortical amyloid-rich differential protein expression (Fig 4F)
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("DESeq2")
  library("ggrepel")
  library("UpSetR")
  library("MAST")
})

# Define output folder
output_folder <- "/path/to/cortical/amyloid/protein/DE/output/folder/"
  
# Load integrated protein Seurat object
s <- readRDS("/path/to/all_samples_03_pro_filtered.rds")

# Load integrated RNA Seurat object to add meta data
temp <- readRDS("/path/to/all_samples_03_rna_updated.rds")

# Define filter operator 
`%notin%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
# Run differential expression 

# Add meta data from RNA object to protein object
sum(row.names(temp@meta.data) != row.names(s@meta.data))
s@meta.data <- cbind(s@meta.data, temp@meta.data[,18:ncol(temp@meta.data)])

# Define amyloid enrichment
s$amyloid_rich <- "not_rich"
s$amyloid_rich[s$plaque_fluo > 153 & s$vessel_fluo == 0] <- "rich"

# Subset for cortical amyloid-rich spots in gray matter
s <- subset(s, amyloid_rich == "rich" & gray_matter != "not_gray")

# Create abbreviated sample ID variable
s@meta.data$sample_de <- str_split_fixed(s@meta.data$sample_id, "\\.", 2)[,2]
s@meta.data$sample_de <- factor(s@meta.data$sample_de)
print(unique(s@meta.data$sample_de))

# Set idents to sample ID
Idents(s) <- "sample_de"

# Combine data across pools and samples
s <- JoinLayers(s)

# Calculate CDR based on isotype-normalized counts 
cdr <- colMeans(GetAssayData(s, assay = "Protein", layer = "counts") > 0)
s@meta.data$cdr <- cdr

# Standardize CDR
s@meta.data$cdr_centered <- as.numeric((s@meta.data$cdr - mean(s@meta.data$cdr)) / sd(s@meta.data$cdr))

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Run FindMarkers with negative binomial test for each comparison (B = Lecanemab, A = CAA control, 1 = FCX, 3 = TCX, 4 = PCX, 9 = HIPP)
comparisons <- c("B1_vs_A1", "B3_vs_A3", "B4_vs_A4", "B9_vs_A9")
for (comparison in comparisons) {
  print(comparison)
  
  # Create output folder for results
  cur_folder <- paste0(output_folder, comparison, "/")
  dir.create(cur_folder, showWarnings = FALSE, recursive = TRUE) 
  
  # Define idents for current comparison
  ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
  ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
  
  # Run FindMarkers, filtering for antibodies expressed in at least 1% of either group (no filtering by LFC)
  results <- FindMarkers(s, assay = "Protein", ident.1 = ident.1, ident.2 = ident.2, test.use = "negbinom", min.pct = 0.01, logfc.threshold = -Inf, 
                         latent.vars = "cdr_centered")
  
  # Adjust p-values and define DE antibodies
  results$BH <- p.adjust(results$p_val, method = "BH")
  results$antibody <- row.names(results)
  results$DE <- "Not DE"
  results$DE[results$avg_log2FC > fc_thresh & results$BH < p_thresh] <- "Upregulated"
  results$DE[results$avg_log2FC < -fc_thresh & results$BH < p_thresh] <- "Downregulated"
  results$DE_antibody <- NA
  results$DE_antibody[results$DE != "Not DE"] <- results$antibody[results$DE != "Not DE"]
  
  print(nrow(results))
  write.csv(results, paste0(cur_folder, "negbinom_cdr_results.csv"))
}

