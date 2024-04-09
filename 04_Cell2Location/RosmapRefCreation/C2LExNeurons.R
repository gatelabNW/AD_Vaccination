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
# Summary: ROSMAP excitatory neuron sampling for spatial reference
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

EN <- readRDS("EN1.rds")

VlnPlot(EN, assay = "RNA", slot = "counts", group.by = "ct", features = "nCount_RNA", raster = FALSE)

celltypes <- unique(EN@meta.data$ct)

EN@meta.data$cellnumber <- sample(c(1:nrow(EN@meta.data)), nrow(EN@meta.data), replace = F)

EN@meta.data$selected <- selected <- 0

table(EN@meta.data$ct)

EN@meta.data$selected[EN@meta.data$ct %in% c("MYO16 EN")] <- 1

celltypes <- celltypes[-3]

for(cell in celltypes){
  
  selected <- sample(EN@meta.data$cellnumber[EN@meta.data$ct == cell], 2000, replace = FALSE)
  
  EN@meta.data$selected[EN@meta.data$cellnumber %in% selected] <- 1
  
}

table(EN@meta.data$ct, EN@meta.data$selected)

newEN <- subset(EN, selected == 1)

VlnPlot(newEN, assay = "RNA", slot = "counts", group.by = "ct", features = "nCount_RNA", raster = FALSE)

saveRDS(newEN, "EN.rds")