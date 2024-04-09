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
# Summary: Microglia differential expression with MAST (Fig 3J, Ext Fig 3H)
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
output_folder <- "/path/to/microglia/per/region/MAST/results/"
  
# Load label transfer Seurat object with final microglia annotations
load("/path/to/s_label_transfer_final")

#-------------------------------------------------------------------------------
# Run differential expression 

# Subset for microglia
s <- subset(s, subtype_clean %in% c("Microglia", "Dividing Microglia"))

# Create abbreviated sample ID variable
s@meta.data$sample_de <- str_split_fixed(s@meta.data$sample_id, "\\.", 3)[,3]
print(unique(s@meta.data$sample_de))

# Set default assay 
DefaultAssay(object = s) <- "SCT"

# Recorrect SCT counts across layers (layers are split by pool and sample ID)
s <- PrepSCTFindMarkers(object = s)

# Calculate CDR using SCT expression
cdr <- colMeans(GetAssayData(s, assay = "SCT", layer = "data") > 0)
s@meta.data$cdr <- cdr

# Standardize CDR and make sample ID factor
s@meta.data$cdr_centered <- as.numeric((s@meta.data$cdr - mean(s@meta.data$cdr)) / sd(s@meta.data$cdr))
s@meta.data$sample_de <- factor(s@meta.data$sample_de)

# Set idents to sample ID
Idents(s) <- "sample_de"

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Extract SCT expression data
expressionmat_full <- GetAssayData(s, assay = "SCT", layer = 'data')
expressionmat_full <- as.matrix(expressionmat_full)

# Run MAST for each comparison (B = Lecanemab, A = CAA control, 1 = FCX, 3 = TCX, 4 = PCX, 9 = HIPP)
comparisons <- c("B1_vs_A1", "B3_vs_A3", "B4_vs_A4", "B9_vs_A9", "B3_vs_B1", "B4_vs_B1", "B9_vs_B1")
for (comparison in comparisons) {
  print(comparison)
  
  # Create output folder for results
  cur_folder <- paste0(output_folder, comparison, "/")
  dir.create(cur_folder, showWarnings = FALSE, recursive = TRUE)
  
  # Define idents for current comparison
  ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
  ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
  
  # Get avg_log2FC and percent expression for all genes
  LFC <- FoldChange(object = s, ident.1 = ident.1, ident.2 = ident.2, fc.name = "avg_log2FC", base = 2) %>% data.frame()
  
  # Filter for genes expressed in at least 1% of either group
  genes_keep <- row.names(LFC)[LFC$pct.1 >= 0.01 | LFC$pct.2 >= 0.01]
  
  # Remove contamination
  if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
    genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)] 
  }
  
  # Filter expression matrix
  expressionmat <- expressionmat_full[row.names(expressionmat_full) %in% genes_keep,]
  
  # Extract cell-level and feature-level meta data for MAST
  cdat <- s@meta.data
  cdat$wellKey <- row.names(cdat)
  fdat <- data.frame(primerid = row.names(expressionmat), row.names = row.names(expressionmat))
  
  # Subset to samples for comparison
  cdat <- cdat[cdat$sample_de %in% c(ident.1, ident.2),]
  expressionmat <- expressionmat[,colnames(expressionmat) %in% row.names(cdat)]
  print(sum(colnames(expressionmat) != row.names(cdat)))
  
  # Create SingleCellAssay object
  sca <- FromMatrix(expressionmat, cdat, fdat, check_sanity = FALSE)
  
  # Set reference level for comparison 
  sample <- factor(colData(sca)$sample_de)
  sample <- relevel(sample, ident.1) # FindMarkers: reference level is "Group1"
  colData(sca)$sample_de <- sample
  
  # Run DE using generalized linear model
  tryCatch(
    {
      zlm_sample <- zlm(~ sample_de + cdr_centered, sca)
      
      lrt_name <- paste0("sample_de", ident.2)
      summary <- summary(zlm_sample, doLRT = lrt_name) # FindMarkers: doLRT = "conditionGroup2"
      
      # Extract hurdle p-values
      summary_data <- summary$datatable %>% data.frame()
      p_val <- summary_data[summary_data[, "component"] == "H", 4]
      genes.return <- summary_data[summary_data[, "component"] == "H", 1]
      
      # Compile results
      results <- data.frame(p_val = p_val, gene = genes.return, row.names = genes.return)
      results$BH <- p.adjust(results$p_val, method = "BH")
      
      # Add LFC data
      LFC_use <- LFC[match(results$gene, row.names(LFC)),]
      print(sum(row.names(LFC_use) != results$gene))
      results$avg_log2FC <- LFC_use$avg_log2FC
      
      # Define DEGs
      results$DE <- "Not DE"
      results$DE[results$avg_log2FC > fc_thresh & results$BH < p_thresh] <- "Upregulated"
      results$DE[results$avg_log2FC < -fc_thresh & results$BH < p_thresh] <- "Downregulated"
      
      # Create variable for volcano plot labels
      results$DE_gene <- NA
      results$DE_gene[results$DE != "Not DE"] <- results$gene[results$DE != "Not DE"]
      
      write.csv(results, paste0(cur_folder, "mast_results.csv"))
    },
    error = function(e) {
      print(paste0("Error in: ", comparison))
    }
  )
}




