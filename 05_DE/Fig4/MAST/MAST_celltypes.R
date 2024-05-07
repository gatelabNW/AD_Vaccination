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
# Summary: Cell type differential expression with MAST (Fig 3I, Ext Fig 3G)
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
output_folder <- "/path/to/celltype/all/region/MAST/results/"
  
# Load label transfer Seurat object with final microglia annotations
load("/path/to/s_label_transfer_final")

# Define filter operator
`%notin%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
# Run differential expression 

# Combine cell type groups for DE
s@meta.data$cell_type_de <- s@meta.data$subtype_clean
s@meta.data$cell_type_de[grep("Microglia", s@meta.data$cell_type_de)] <- "Microglia"

# Filter for cells with subtype annotation 
s <- subset(s, subtype_clean != "Unknown")

# Create merged group variable ("A" = CAA control samples, "B" = Lecanemab samples)
s@meta.data$sample_de <- str_split_fixed(s@meta.data$sample_id, "\\.", 3)[,3]
s@meta.data$group_de <- NA
s@meta.data$group_de[str_detect(s@meta.data$sample_de, "B")] <- "B"
s@meta.data$group_de[str_detect(s@meta.data$sample_de, "A")] <- "A"
print(unique(s@meta.data$group_de))

# Make group and sample factors
s@meta.data$group_de <- factor(s@meta.data$group_de)
s@meta.data$sample_de <- factor(s@meta.data$sample_de)

# Set idents to group ID 
Idents(s) <- "group_de"

# Set default assay
DefaultAssay(object = s) <- "SCT"

# Create output folder for results
comparison <- "B_vs_A"
result_folder <- paste0(output_folder, comparison, "/")
dir.create(result_folder, showWarnings = FALSE, recursive = TRUE)

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) 

# Run MAST for each cell type
for (type in unique(s@meta.data$cell_type_de)) {
  print(type)
  
  # Remove samples with less than 10 cells 
  meta <- s@meta.data[s@meta.data$cell_type_de == type,]
  n_cells <- table(meta$sample_de)
  samples_keep <- names(n_cells[n_cells >= 10])
  if (length(samples_keep) != length(names(n_cells))) {
    print("Samples removed: ")
    print(names(n_cells[n_cells < 10]))
  }
  
  # Subset for ROI
  s_type <- subset(s, cell_type_de == type & sample_de %in% samples_keep)
  
  # Recorrect SCT counts across layers (layers are split by pool and sample ID)
  s_type <- PrepSCTFindMarkers(object = s_type)
  
  # Calculate CDR using SCT expression
  cdr <- colMeans(GetAssayData(s_type, assay = "SCT", layer = "data") > 0)
  s_type@meta.data$cdr <- cdr
  
  # Standardize CDR within ROI
  s_type@meta.data$cdr_centered <- as.numeric((s_type@meta.data$cdr - mean(s_type@meta.data$cdr)) / sd(s_type@meta.data$cdr))
  
  # Extract SCT expression data
  expressionmat_full <- GetAssayData(s_type, assay = "SCT", layer = 'data')
  expressionmat_full <- as.matrix(expressionmat_full)
  
  # Define idents for current comparison
  ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
  ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
  
  # Confirm cell type is present in both groups
  if (ident.1 %notin% unique(Idents(s_type)) | ident.2 %notin% unique(Idents(s_type))) {
    print("Both conditions not present")
    next
  }
  
  # Get avg_log2FC and percent expression for all genes
  LFC <- FoldChange(object = s_type, ident.1 = ident.1, ident.2 = ident.2, fc.name = "avg_log2FC", base = 2) %>% data.frame()
  
  # Filter for genes expressed in a minimum percent of either group (threshold adjusted for certain analyses - MAST cannot fit some models due to low expression)
  if (type %in% c("FB", "Oligodendrocyte", "Pericyte", "OPC")) {
    genes_keep <- row.names(LFC)[LFC$pct.1 >= 0.15 | LFC$pct.2 >= 0.15]
  } else {
    genes_keep <- row.names(LFC)[LFC$pct.1 >= 0.1 | LFC$pct.2 >= 0.1]
  }
  
  # Remove contamination
  if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
    genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)] 
  }
  
  # Filter expression matrix
  expressionmat <- expressionmat_full[row.names(expressionmat_full) %in% genes_keep,]
  
  # Extract cell-level and feature-level meta data for MAST
  cdat <- s_type@meta.data
  cdat$wellKey <- row.names(cdat)
  fdat <- data.frame(primerid = row.names(expressionmat), row.names = row.names(expressionmat))
  
  # Make sure columns of expression data match rows of meta data
  print(sum(colnames(expressionmat) != row.names(cdat)))
  
  # Create SingleCellAssay object
  sca <- FromMatrix(expressionmat, cdat, fdat, check_sanity = FALSE)
  
  # Set reference level for comparison 
  group <- factor(colData(sca)$group_de)
  group <- relevel(group, ident.1) # FindMarkers: reference level is "Group1"
  colData(sca)$group_de <- group
  
  # Run DE using generalized linear mixed model
  tryCatch(
    {
      # MAST does not allow ebayes (empirical Bayes estimation to regularize variance) for glmer
      zlm_group <- zlm(~ group_de + cdr_centered + (1 | sample_de), sca, ebayes = FALSE, method = "glmer")
      
      lrt_name <- paste0("group_de", ident.2)
      summary <- summary(zlm_group, doLRT = lrt_name) # FindMarkers: doLRT = "conditionGroup2"
      
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
      
      type_name <- str_replace_all(type, " ", "_")
      type_name <- str_replace_all(type_name, "/", "_")
      write.csv(results, paste0(result_folder, type_name, "_mast_results.csv"))
    },
    error = function(e) {
      print(paste0("Error in: ", type))
    }
  )
}