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
# Summary: Pathway analysis of DESeq2 DEGs using gseKEGG 
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(pathview)
  library(ggplot2)
  library(dplyr)
})

# Define output folder
output_folder <- "/path/to/general/cohort1/DE/output/folder/"

#-------------------------------------------------------------------------------
# Run gseKEGG for gray matter microglia: nAD vs. NNC, lim vs. nAD, ext vs. nAD 

# Load data
nAD_NNC <- read.csv(paste0(output_folder, "cell_types_gray/nAD_vs_NNC/deseq/results/Microglia_T5.csv"), row.names = 1)
lim_nAD <- read.csv(paste0(output_folder, "cell_types_gray/lim_vs_nAD/deseq/results/Microglia_T5.csv"), row.names = 1)
ext_nAD <- read.csv(paste0(output_folder, "cell_types_gray/ext_vs_nAD/deseq/results/Microglia_T5.csv"), row.names = 1)

# Get gene symbols
nAD_NNC_genes <- nAD_NNC$gene
lim_nAD_genes <- lim_nAD$gene
ext_nAD_genes <- ext_nAD$gene

nAD_NNC_entrez <- bitr(nAD_NNC_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
lim_nAD_entrez <- bitr(lim_nAD_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
ext_nAD_entrez <- bitr(ext_nAD_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# Ensure that there is only one entrez ID per gene symbol
nAD_NNC_entrez <- nAD_NNC_entrez[match(unique(nAD_NNC_entrez$SYMBOL), nAD_NNC_entrez$SYMBOL),]
lim_nAD_entrez <- lim_nAD_entrez[match(unique(lim_nAD_entrez$SYMBOL), lim_nAD_entrez$SYMBOL),]
ext_nAD_entrez <- ext_nAD_entrez[match(unique(ext_nAD_entrez$SYMBOL), ext_nAD_entrez$SYMBOL),]

# Filter for genes mapped to entrez ID
nAD_NNC <- nAD_NNC[nAD_NNC$gene %in% nAD_NNC_entrez$SYMBOL,]
lim_nAD <- lim_nAD[lim_nAD$gene %in% lim_nAD_entrez$SYMBOL,]
ext_nAD <- ext_nAD[ext_nAD$gene %in% ext_nAD_entrez$SYMBOL,]

# Confirm order matches before adding entrez IDs
nAD_NNC$X <- nAD_NNC$gene
lim_nAD$X <- lim_nAD$gene
ext_nAD$X <- ext_nAD$gene

sum(nAD_NNC$X != nAD_NNC_entrez$SYMBOL)
sum(lim_nAD$X != lim_nAD_entrez$SYMBOL)
sum(ext_nAD$X != ext_nAD_entrez$SYMBOL)

nAD_NNC$Y <- nAD_NNC_entrez$ENTREZID
lim_nAD$Y <- lim_nAD_entrez$ENTREZID
ext_nAD$Y <- ext_nAD_entrez$ENTREZID

# Calculate PFC 
nAD_NNC$PFC <- -log10(nAD_NNC$padj)*nAD_NNC$log2FoldChange
lim_nAD$PFC <- -log10(lim_nAD$padj)*lim_nAD$log2FoldChange
ext_nAD$PFC <- -log10(ext_nAD$padj)*ext_nAD$log2FoldChange

# Create a vector of the gene universe
nAD_NNC_kegg_genes <- nAD_NNC$PFC
lim_nAD_kegg_genes <- lim_nAD$PFC
ext_nAD_kegg_genes <- ext_nAD$PFC

names(nAD_NNC_kegg_genes) <- nAD_NNC$Y
names(lim_nAD_kegg_genes) <- lim_nAD$Y
names(ext_nAD_kegg_genes) <- ext_nAD$Y

nAD_NNC_kegg_genes <- sort(nAD_NNC_kegg_genes, decreasing = TRUE)
lim_nAD_kegg_genes <- sort(lim_nAD_kegg_genes, decreasing = TRUE)
ext_nAD_kegg_genes <- sort(ext_nAD_kegg_genes, decreasing = TRUE)

# Set organism to human
kegg_organism <- 'hsa'

# Generate test results with multiple seeds for each comparison to ensure results are consistent
for (seed in c(50, 100, 150, 200)) {
  set.seed(seed)
  # results <- gseKEGG(geneList = nAD_NNC_kegg_genes, organism = kegg_organism, pvalueCutoff = 0.05, seed = TRUE,
  #                    pAdjustMethod = "BH", keyType = "ncbi-geneid")
  # results <- gseKEGG(geneList = lim_nAD_kegg_genes, organism = kegg_organism, pvalueCutoff = 0.05, seed = TRUE,
  #                    pAdjustMethod = "BH", keyType = "ncbi-geneid")
  results <- gseKEGG(geneList = ext_nAD_kegg_genes, organism = kegg_organism, pvalueCutoff = 0.05, seed = TRUE,
                     pAdjustMethod = "BH", keyType = "ncbi-geneid")
  print(results@result$Description)
}

# Set seed for final results to ensure reproducibility
set.seed(100)

# Generate final results for each group
nAD_NNC_kegg_results <- gseKEGG(geneList = nAD_NNC_kegg_genes, organism = kegg_organism, pvalueCutoff = 0.05, seed = TRUE,
                                pAdjustMethod = "BH", keyType = "ncbi-geneid")

lim_nAD_kegg_results <- gseKEGG(geneList = lim_nAD_kegg_genes, organism = kegg_organism, pvalueCutoff = 0.05, seed = TRUE,
                                pAdjustMethod = "BH", keyType = "ncbi-geneid")

ext_nAD_kegg_results <- gseKEGG(geneList = ext_nAD_kegg_genes, organism = kegg_organism, pvalueCutoff = 0.05, seed = TRUE,
                                pAdjustMethod = "BH", keyType = "ncbi-geneid")

# Save results
saveRDS(nAD_NNC_kegg_results, paste0(output_folder, "custom_output/deseq_unfiltered/gseKEGG_all_genes/microglia_gray/nAD_NNC_adj.rds"))
saveRDS(lim_nAD_kegg_results, paste0(output_folder, "custom_output/deseq_unfiltered/gseKEGG_all_genes/microglia_gray/lim_nAD_adj.rds"))
saveRDS(ext_nAD_kegg_results, paste0(output_folder, "custom_output/deseq_unfiltered/gseKEGG_all_genes/microglia_gray/ext_nAD_adj.rds"))


