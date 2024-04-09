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
# Summary: ROSMAP oligodendrocyte precursor cell sampling for spatial reference
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

OPC <- readRDS("OPC_annotation.rds")

VlnPlot(OPC, assay = "RNA", slot = "counts", group.by = "ct", features = "nCount_RNA")

table(OPC@meta.data$ct)

OPC@meta.data$cellnumber <- sample(c(1:nrow(OPC@meta.data)), nrow(OPC@meta.data), replace = F)

OPC@meta.data$selected <- selected <- 0

selected <- sample(OPC@meta.data$cellnumber, 2000, replace = FALSE)

OPC@meta.data$selected[OPC@meta.data$cellnumber %in% selected] <- 1

table(OPC@meta.data$ct, OPC@meta.data$selected)

newOPC <- subset(OPC, selected == 1)

newOPC$ct <- "OPC"

saveRDS(newOPC, "OPC.rds")