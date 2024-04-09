# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AN1792 Project                              -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Written by: Thomas Watson
# Summary: Label spots for cell-type enrichment using C2L predicted abundances
#
#-------------------------------------------------------------------------------

suppressMessages({
  library(plyr)
  library(data.table)
  library(tidyverse) 
  library(Seurat)
  library(rlist)
  library(doMC)
  library(UpSetR)
  library(DElegate)
  library(DESeq2)
  library(edgeR)
  library(sparseMatrixStats)
  library(SeuratData)
  library(ComplexHeatmap)
  library(ggrepel)
  library(rsvd)
  library(plotly)
  library(ggthemes)
  library(ggplot2)
  library(cowplot)
  library(ggsci)
  library(ggdark)
  library(viridis)
  library(car)
  library(e1071)
  library(parallel)
  library(knitr)
  library(rlist)
  library(SCP)
  library(DescTools)
  library(ppcor)
  library(fastDummies)
})

options(future.globals.maxSize=1048576000000)
load("helperfunctions.RData")
setwd("/folder/containing/C2L/output")

# Load seurat object used for C2L
s <- readRDS("seuratobject.rds")

# get filenames of per-sample metas from 11/21 C2L run
fnames <- list.files("C2L/Output")[-27]

# initialize list for metas
metas <- list()

# Get names in right format
IDs <- gsub("_.*", "", fnames)

# Populate metas with C2L output
for(i in 1:length(fnames)){
  
  metas[[IDs[i]]] <- read.csv(paste0(fnames[[i]]))
  
}

# Bind results list
newmeta <- list.rbind(metas)

# 
s@meta.data <- s@meta.data[s@meta.data$sample_barcode %in% newmeta$sample_barcode,]

newmeta <- newmeta[newmeta$sample_barcode %in% s@meta.data$sample_barcode,]

# some preds are NA like oligos in meninges etc - make them 0
newmeta[is.na(newmeta)] <- 0

# order by barcode
newmeta <- newmeta[order(newmeta$sample_barcode),]

# order seurat object meta by barcode
s@meta.data <- s@meta.data[order(s@meta.data$sample_barcode),]

# make sure it's safe to attach C2L preds
all.equal(newmeta$sample_barcode, s@meta.data$sample_barcode)

# remove old preds
s@meta.data <- s@meta.data[,c(1:33)]

# attach C2L preds - we already have barcode
s@meta.data <- cbind(s@meta.data, newmeta[,25:42])

# Set manual annotation as idents for plots
Idents(s) = s@meta.data$manual_annotation

# Loop over samples
for(sample in unique(s@meta.data$sample_id)){
  
  # manual layer, sample_barcode, and C2L preds
  T5 <- metas[[sample]][,c(11, 13, 25:42)]
  
  # so names are different when we cbind preds to meta data
  colnames(T5) <- paste0(colnames(T5), "_T5")
  
  # G1
  
  T5$Pericyte_T5 <- ifelse(T5$Pericyte_T5 > quantile(T5$Pericyte_T5 , probs = 0.95), 1, 0)
  
  T5$Peripheral.Immune.Cell_T5 <- ifelse(T5$Peripheral.Immune.Cell_T5 > quantile(T5$Peripheral.Immune.Cell_T5 , probs = 0.95), 1, 0)
  
  T5$SMC_T5 <- ifelse(T5$SMC_T5 > quantile(T5$SMC_T5 , probs = 0.95), 1, 0)
  
  T5$Endothelial.Cell_T5 <- ifelse(T5$Endothelial.Cell_T5 > quantile(T5$Endothelial.Cell_T5 , probs = 0.95), 1, 0)
  
  T5$FB_T5 <- ifelse(T5$FB_T5 > quantile(T5$FB_T5 , probs = 0.95), 1, 0)
  
  # G2
  
  T5$Interneuron_T5[T5$manual_annotation_T5 %in% c("meninges", "white")] <- 0
  
  T5$Interneuron_T5 <- ifelse(T5$Interneuron_T5 > quantile(T5$Interneuron_T5 , probs = 0.95), 1, 0)
  
  T5$MYO16.EN_T5[T5$manual_annotation_T5 %in% c("meninges", "white")] <- 0
  
  T5$MYO16.EN_T5 <- ifelse(T5$MYO16.EN_T5 > quantile(T5$MYO16.EN_T5 , probs = 0.95), 1, 0)
  
  # G3
  
  T5$Oligodendrocyte_T5[T5$manual_annotation_T5 %in% c("meninges")] <- 0
  
  T5$Oligodendrocyte_T5 <- ifelse(T5$Oligodendrocyte_T5 > quantile(T5$Oligodendrocyte_T5 , probs = 0.95), 1, 0)
  
  T5$Astrocyte_T5[T5$manual_annotation_T5 %in% c("meninges")] <- 0
  
  T5$Astrocyte_T5 <- ifelse(T5$Astrocyte_T5 > quantile(T5$Astrocyte_T5 , probs = 0.95), 1, 0)
  
  T5$OPC_T5[T5$manual_annotation_T5 %in% c("meninges")] <- 0
  
  T5$OPC_T5 <- ifelse(T5$OPC_T5 > quantile(T5$OPC_T5 , probs = 0.95), 1, 0)
  
  # G4
  
  glia <- data.frame("sample_barcode" = T5$sample_barcode_T5, "layer" = T5$manual_annotation_T5, "glia" = T5$Microglia_T5)
  
  glia$glia[glia$layer == "meninges"] <- 0
  
  glia$glia[glia$layer == "white"] <- ifelse(glia$glia[glia$layer == "white"] > quantile(glia$glia[glia$layer == "white"], probs = 0.95), 1, 0)
  
  glia$glia[glia$layer != "white"] <- ifelse(glia$glia[glia$layer != "white"] > quantile(glia$glia[glia$layer != "white"], probs = 0.95), 1, 0)
  
  print(all.equal(glia$sample_barcode, T5$sample_barcode_T5))
  
  T5$Microglia_T5 <- glia$glia
  
  # G5
  
  T5$L2.3.EN_T5[T5$manual_annotation_T5 %!in% c("gray-l2", "gray-l3")] <- 0
  
  T5$L2.3.EN_T5 <- ifelse(T5$L2.3.EN_T5 > quantile(T5$L2.3.EN_T5 , probs = 0.95), 1, 0)
  
  T5$L4.EN_T5[T5$manual_annotation_T5 %!in% c("gray-l4")] <- 0
  
  T5$L4.EN_T5 <- ifelse(T5$L4.EN_T5 > quantile(T5$L4.EN_T5 , probs = 0.95), 1, 0)
  
  T5$L4.5.EN_T5[T5$manual_annotation_T5 %!in% c("gray-l4", "gray-l56")] <- 0
  
  T5$L4.5.EN_T5 <- ifelse(T5$L4.5.EN_T5 > quantile(T5$L4.5.EN_T5 , probs = 0.95), 1, 0)
  
  T5$L5.EN_T5[T5$manual_annotation_T5 %!in% c("gray-l56")] <- 0
  
  T5$L5.EN_T5 <- ifelse(T5$L5.EN_T5 > quantile(T5$L5.EN_T5 , probs = 0.95), 1, 0)
  
  T5$L5.6.EN_T5[T5$manual_annotation_T5 %!in% c("gray-l56")] <- 0
  
  T5$L5.6.EN_T5 <- ifelse(T5$L5.6.EN_T5 > quantile(T5$L5.6.EN_T5 , probs = 0.95), 1, 0)
  
  T5$L5.6.CCa.EN_T5[T5$manual_annotation_T5 %!in% c("gray-l56")] <- 0
  
  T5$L5.6.CCa.EN_T5 <- ifelse(T5$L5.6.CCa.EN_T5 > quantile(T5$L5.6.CCa.EN_T5 , probs = 0.95), 1, 0)
  
  T5$L5.6.CCb.EN_T5[T5$manual_annotation_T5 %!in% c("gray-l56")] <- 0
  
  T5$L5.6.CCb.EN_T5 <- ifelse(T5$L5.6.CCb.EN_T5 > quantile(T5$L5.6.CCb.EN_T5 , probs = 0.95), 1, 0)
  
  print(all.equal(metas[[sample]]$sample_barcode, T5$sample_barcode_T5))
  
  metas[[sample]] <- cbind(metas[[sample]], T5[,3:20])
  
}

newmetas <- list.rbind(metas)

newmetas <- newmetas[newmetas$sample_barcode %in% s@meta.data$sample_barcode,]

newmetas <- newmetas[order(newmetas$sample_barcode),]

all.equal(s@meta.data$sample_barcode, newmetas$sample_barcode)

s@meta.data <- cbind(s@meta.data, newmetas[,45:62])

snames <- unique(s@meta.data$sample_id)



cts<- colnames(s@meta.data[,34:51])


# Looping over sample ID 
for(sample in unique(s@meta.data$sample_id)){
  
  # Open PDF to store 18 celltype plots per sample
  pdf(file = paste0("/C2L/Output/Plots/", sample, "_C2Lpreds.pdf"), width = 8, height = 8)
  
  # subset to sample
  s2 <- subset(s, sample_id == sample)
  
  # now loop over scaled celltypes
  for(celltype in cts){
    
    # Plot scaled celltype preds, split by manual layer
    p <- VlnPlot(s2, features = celltype, raster = FALSE) + facet_grid(.~"manual_annotation") +
      ggtitle(paste0(celltype)) + xlab(sample)
    
    p$data$ident <- factor(x = p$data$ident, levels = c("meninges", "gray-l1", "gray-l2", "gray-l3", "gray-l4", "gray-l56", "white", ""))
    
    # Some spots don't have manual annotation? Removing them
    p$data <- p$data[p$data$ident != "",]
    
    # Print this celltype's plot in the sample's pdf
    print(p)
    
  }
  
  # Finish writing PDF so we can start a new one for the next sample
  dev.off()
  
}