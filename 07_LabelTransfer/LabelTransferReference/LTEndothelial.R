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
# Summary: ROSMAP endothelial cell sampling for label transfer reference
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

endothelial <- readRDS("endothelial_cell_annotation.rds")

VlnPlot(endothelial, assay = "RNA", slot = "counts", group.by = "ct", features = "nCount_RNA")

endothelial <- subset(endothelial, ct %!in% c("EC-like pericyte", "EC-like SMC", "Immune-EC doublet"))

table(endothelial@meta.data$ct)

endothelial@meta.data$cellnumber <- sample(c(1:nrow(endothelial@meta.data)), nrow(endothelial@meta.data), replace = F)

endothelial@meta.data$selected <- selected <- 0

selected <- sample(endothelial@meta.data$cellnumber, 5000, replace = FALSE)

endothelial@meta.data$selected[endothelial@meta.data$cellnumber %in% selected] <- 1

table(endothelial@meta.data$ct, endothelial@meta.data$selected)

newendothelial <- subset(endothelial, selected == 1)

newendothelial$ct <- "Endothelial Cell"

saveRDS(newendothelial, "endothelial0129.rds")
