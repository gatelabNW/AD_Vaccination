# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AN1792 Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Written by: Thomas Watson
# Summary: ROSMAP stromal sampling for label transfer reference
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
  library(aricode)
  library(knitr)
  library(rlist)
  library(SCP)
})

options(future.globals.maxSize=1048576000000)
load("helperfunctions.RData")
setwd("/objects/")
set.seed(256)

stromal <- readRDS("stromal_cell_annotation.rds")

stromal2 <- readRDS("counts/stromal_cell_annotation.rds")

newcounts <- GetAssayData(stromal2, layer = "counts", assay = "raw")

stromal@assays$RNA$counts <- newcounts

VlnPlot(stromal, assay = "RNA", slot = "counts", group.by = "ct", features = "nCount_RNA")

stromal <- subset(stromal, ct != "AST-FB doublet")

celltypes <- unique(stromal@meta.data$ct)

stromal@meta.data$cellnumber <- sample(c(1:nrow(stromal@meta.data)), nrow(stromal@meta.data), replace = F)

stromal@meta.data$selected <- selected <- 0

table(stromal@meta.data$ct)

stromal@meta.data$selected[stromal@meta.data$ct == "SMC"] <- 1

celltypes <- celltypes[1:2]

for(cell in celltypes){
  
  selected <- sample(stromal@meta.data$cellnumber[stromal@meta.data$ct == cell], 5000, replace = FALSE)
  
  stromal@meta.data$selected[stromal@meta.data$cellnumber %in% selected] <- 1
  
}

table(stromal@meta.data$ct, stromal@meta.data$selected)

newStromal <- subset(stromal, selected == 1)

View(newStromal@meta.data)

saveRDS(newStromal, "stromal0129.rds")