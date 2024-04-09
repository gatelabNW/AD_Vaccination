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
# Summary: ROSMAP immune cell sampling for spatial reference
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

immune <- readRDS("immune_cell_annotation.rds")

VlnPlot(immune, assay = "RNA", slot = "counts", group.by = "ct", features = "nCount_RNA")

table(immune@meta.data$ct)

immune@meta.data$ct[immune@meta.data$ct %!in% c("microglia")] <- "Peripheral Immune Cell"

immune@meta.data$ct[immune@meta.data$ct %in% c("microglia")] <- "Microglia"

celltypes <- unique(immune@meta.data$ct)

immune@meta.data$cellnumber <- sample(c(1:nrow(immune@meta.data)), nrow(immune@meta.data), replace = F)

immune@meta.data$selected <- selected <- 0

for(cell in celltypes){
  
  selected <- sample(immune@meta.data$cellnumber[immune@meta.data$ct == cell], 2000, replace = FALSE)
  
  immune@meta.data$selected[immune@meta.data$cellnumber %in% selected] <- 1
  
}

table(immune@meta.data$ct, immune@meta.data$selected)

newimmune <- subset(immune, selected == 1)

table(newimmune$ct)

saveRDS(newimmune, "immune.rds")

ImmuneMNG <- subset(newimmune, ct != "Microglia")

table(ImmuneMNG$ct)

saveRDS(ImmuneMNG, "MNGImmune.rds")