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
# Summary: ROSMAP immune cell sampling for label transfer reference
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

interneuron <- readRDS("interneuron_annotation.rds")

VlnPlot(interneuron, assay = "RNA", slot = "counts", group.by = "ct", features = "nFeature_RNA")

table(interneuron@meta.data$ct)

interneuron@meta.data$cellnumber <- sample(c(1:nrow(interneuron@meta.data)), nrow(interneuron@meta.data), replace = F)

interneuron@meta.data$selected <- selected <- 0

selected <- sample(interneuron@meta.data$cellnumber, 5000, replace = FALSE)

interneuron@meta.data$selected[interneuron@meta.data$cellnumber %in% selected] <- 1

table(interneuron@meta.data$ct, interneuron@meta.data$selected)

newinterneuron <- subset(interneuron, selected == 1)

newinterneuron$ct <- "Interneuron"

saveRDS(newinterneuron, "interneuron0129.rds")