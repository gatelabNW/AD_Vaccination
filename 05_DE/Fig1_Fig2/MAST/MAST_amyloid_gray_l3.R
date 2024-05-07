# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                         AN1792 Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 03-29-2024
# Written by: Anne Forsyth
# Summary: Differential expression of amyloid proximity groups in gray matter layer 3 (Fig 2I)
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
output_folder <- "/path/to/amyloid/proximity/MAST/results/"

# Load integrated Seurat object
s <- readRDS("/path/to/all_samples_03.rds")

# Define filter operator
`%notin%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
# Run differential expression 

# Remove spots to exclude from amyloid analysis
s <- subset(s, amyloid_filter == "include")

# Subset for gray matter layer 3
s <- subset(s, manual_annotation == "gray-l3")

# Define broad amyloid proximity groups
s@meta.data <- s@meta.data %>% mutate(amyloid_de = case_when(amyloid_distance == "amyloid_0" ~ "plaque_far", 
                                                             amyloid_distance %in% c("amyloid_1", "amyloid_2", "amyloid_3") ~ "plaque_nearby", 
                                                             amyloid_distance %in% c("amyloid_4", "amyloid_5") ~ "plaque_rich"))

# Identify samples with at least 10 spots for each amyloid group, using all NNC spots
far <- s@meta.data[s@meta.data$amyloid_de == "plaque_far",]
near <- s@meta.data[s@meta.data$amyloid_de == "plaque_nearby",]
rich <- s@meta.data[s@meta.data$amyloid_de == "plaque_rich",]

far_samples <- table(far$sample_id[far$condition != "NNC"])
near_samples <- table(near$sample_id[near$condition != "NNC"]) 
rich_samples <- table(rich$sample_id[rich$condition != "NNC"])
nnc_samples <- table(s@meta.data$sample_id[s@meta.data$condition == "NNC"]) 

# List of samples to keep for each analysis 
samples_keep <- list(names(far_samples[far_samples >= 10]), names(near_samples[near_samples >= 10]),
                     names(rich_samples[rich_samples >= 10]), names(nnc_samples[nnc_samples >= 10]))
names(samples_keep) <- c("plaque_far", "plaque_nearby", "plaque_rich", "NNC")

# Set default assay
DefaultAssay(object = s) <- "SCT"

# Define DE thresholds
p_thresh <- 0.05
fc_thresh <- log2(1.5) # 50% fold change

# Run MAST for "nearby" and "rich" groups
for (cluster in c("plaque_nearby", "plaque_rich")) {
  print(cluster)
  
  # Subset object to current cluster and recorrect SCT data - keep all NNC spots like we did for DESeq2 (to keep ROI the same)
  clust_s <- subset(s, (amyloid_de == cluster | condition == "NNC"))
  
  # Keep samples with at least 10 spots in ROI
  clust_s <- subset(clust_s, sample_id %in% c(samples_keep[[cluster]], samples_keep[['NNC']]))
  print(table(clust_s@meta.data$sample_id))
  
  # Recorrect SCT counts across layers (layers are split by sample ID)
  for (name in names(clust_s@assays$SCT@SCTModel.list)) {
    clust_s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
  }
  clust_s <- PrepSCTFindMarkers(object = clust_s)
  
  # Calculate CDR using SCT expression
  cdr <- colMeans(GetAssayData(clust_s, assay = "SCT", layer = "data") > 0)
  clust_s@meta.data$cdr <- cdr
  
  # Standardize continuous covariates and make categorical covariates factors
  clust_s@meta.data$age_centered <- as.numeric((clust_s@meta.data$age - mean(clust_s@meta.data$age)) / sd(clust_s@meta.data$age))
  clust_s@meta.data$cdr_centered <- as.numeric((clust_s@meta.data$cdr - mean(clust_s@meta.data$cdr)) / sd(clust_s@meta.data$cdr))
  clust_s@meta.data$gDNA_centered <- as.numeric((clust_s@meta.data$gDNA_percent - mean(clust_s@meta.data$gDNA_percent)) / sd(clust_s@meta.data$gDNA_percent))
  clust_s@meta.data$sex <- factor(clust_s@meta.data$sex)
  clust_s@meta.data$sample_id <- factor(clust_s@meta.data$sample_id)
  
  # Extract SCT expression data
  expressionmat_full <- GetAssayData(clust_s, assay = "SCT", layer = 'data')
  expressionmat_full <- as.matrix(expressionmat_full)
  
  # Run MAST for each comparison
  for (comparison in c("ext_vs_nAD", "lim_vs_nAD")) {
    
    # Define idents for current comparison
    ident.1 <- str_split_fixed(comparison, "_vs_", 2)[1]
    ident.2 <- str_split_fixed(comparison, "_vs_", 2)[2]
    
    # If comparing iAD-Extensive or iAD-Limited, set idents to condition_clearance
    if ((ident.1 %in% c("ext", "lim")) | (ident.2 %in% c("ext", "lim"))) {
      Idents(clust_s) <- "condition_clearance"
    } else {
      Idents(clust_s) <- "condition"
    }
    
    # Skip comparison if both conditions are not present for the current analysis
    if ((ident.1 %notin% unique(Idents(clust_s))) | (ident.2 %notin% unique(Idents(clust_s)))) {
      print(paste0("Both conditions not present for ", comparison, ": ", cluster))
      next
    }
    
    # Get avg_log2FC and percent expression for all genes
    LFC <- FoldChange(object = clust_s, ident.1 = ident.1, ident.2 = ident.2, fc.name = "avg_log2FC", base = 2) %>% data.frame()
    
    # Filter for genes expressed in a minimum percent of either group (threshold adjusted for certain analyses - MAST cannot fit some models due to low expression)
    if (cluster == "plaque_nearby") {
      genes_keep <- row.names(LFC)[LFC$pct.1 >= 0.1 | LFC$pct.2 >= 0.1]
    } else {
      genes_keep <- row.names(LFC)[LFC$pct.1 >= 0.05 | LFC$pct.2 >= 0.05]
    }
    
    # Remove contamination
    if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
      genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)] 
    }
    
    # Filter expression matrix
    expressionmat <- expressionmat_full[row.names(expressionmat_full) %in% genes_keep,]
    
    # Extract cell-level and feature-level meta data for MAST
    cdat <- clust_s@meta.data
    cdat$wellKey <- row.names(cdat)
    
    # Override condition variable if comparing iAD-Extensive or iAD-Limited
    if ((ident.1 %in% c("ext", "lim")) | (ident.2 %in% c("ext", "lim"))) {
      cdat$condition <- cdat$condition_clearance
    } 
    
    fdat <- data.frame(primerid = row.names(expressionmat), row.names = row.names(expressionmat))
    
    # Subset to conditions for comparison
    cdat <- cdat[cdat$condition %in% c(ident.1, ident.2),]
    expressionmat <- expressionmat[,colnames(expressionmat) %in% row.names(cdat)]
    
    # Create SingleCellAssay object
    sca <- FromMatrix(expressionmat, cdat, fdat, check_sanity = FALSE)
    
    # Set reference level for condition
    cond <- factor(colData(sca)$condition)
    cond <- relevel(cond, ident.1) # FindMarkers: reference level is "Group1"
    colData(sca)$condition <- cond
    
    # Run DE using generalized linear mixed model
    tryCatch(
      {
        # MAST does not allow ebayes (empirical Bayes estimation to regularize variance) for glmer
        zlm_condition <- zlm(~ condition + sex + age_centered + cdr_centered + gDNA_centered + (1 | sample_id), sca, 
                             ebayes = FALSE, method = "glmer")
        
        lrt_name <- paste0("condition", ident.2)
        summary_condition <- summary(zlm_condition, doLRT = lrt_name) # FindMarkers: doLRT = "conditionGroup2"
        
        # Extract hurdle p-values
        summary_data <- summary_condition$datatable %>% data.frame()
        p_val <- summary_data[summary_data[, "component"] == "H", 4]
        genes.return <- summary_data[summary_data[, "component"] == "H", 1]
        
        # Compile results
        results <- data.frame(p_val = p_val, gene = genes.return, row.names = genes.return)
        results$BH <- p.adjust(results$p_val, method = "BH")
        
        # Add LFC data
        LFC_use <- LFC[match(results$gene, row.names(LFC)),]
        results$avg_log2FC <- LFC_use$avg_log2FC
        
        # Define DEGs
        results$DE <- "Not DE"
        results$DE[results$avg_log2FC > fc_thresh & results$BH < p_thresh] <- "Upregulated"
        results$DE[results$avg_log2FC < -fc_thresh & results$BH < p_thresh] <- "Downregulated"
        
        # Create variable for volcano plot labels
        results$DE_gene <- NA
        results$DE_gene[results$DE != "Not DE"] <- results$gene[results$DE != "Not DE"]
        
        write.csv(results, paste0(output_folder, comparison, "/mast/results/", cluster, ".csv"))
      },
      error = function(e) {
        print(paste0("Error in ", cluster, ": ", comparison))
      }
    )
  }
}










