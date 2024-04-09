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
# Summary: ROSMAP astrocyte sampling for label transfer reference
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

astrocyte <- readRDS("astrocyte_annotation.rds")

VlnPlot(astrocyte, assay = "RNA", slot = "counts", group.by = "ct", features = "nFeature_RNA")

table(astrocyte@meta.data$ct)

celltypes <- unique(astrocyte@meta.data$ct)

astrocyte@meta.data$cellnumber <- sample(c(1:nrow(astrocyte@meta.data)), nrow(astrocyte@meta.data), replace = F)

astrocyte@meta.data$selected <- selected <- 0

selected <- sample(astrocyte@meta.data$cellnumber, 5000, replace = FALSE)

astrocyte@meta.data$selected[astrocyte@meta.data$cellnumber %in% selected] <- 1

table(astrocyte@meta.data$ct, astrocyte@meta.data$selected)

newastrocyte <- subset(astrocyte, selected == 1)

newastrocyte$ct <- "Astrocyte"

saveRDS(newastrocyte, "astrocyte0129.rds")