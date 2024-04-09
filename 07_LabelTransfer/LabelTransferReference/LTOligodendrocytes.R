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
# Summary: ROSMAP oligodendrocyte sampling for label transfer reference
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

oligodendrocyte <- readRDS("oligodendrocyte_annotation.rds")

VlnPlot(oligodendrocyte, assay = "RNA", slot = "counts", group.by = "ct", features = "nCount_RNA")

table(oligodendrocyte@meta.data$ct)

oligodendrocyte@meta.data$cellnumber <- sample(c(1:nrow(oligodendrocyte@meta.data)), nrow(oligodendrocyte@meta.data), replace = F)

oligodendrocyte@meta.data$selected <- selected <- 0

selected <- sample(oligodendrocyte@meta.data$cellnumber, 5000, replace = FALSE)

oligodendrocyte@meta.data$selected[oligodendrocyte@meta.data$cellnumber %in% selected] <- 1

table(oligodendrocyte@meta.data$ct, oligodendrocyte@meta.data$selected)

newoligodendrocyte <- subset(oligodendrocyte, selected == 1)

newoligodendrocyte$ct <- "Oligodendrocyte"

saveRDS(newoligodendrocyte, "oligodendrocyte0129.rds")