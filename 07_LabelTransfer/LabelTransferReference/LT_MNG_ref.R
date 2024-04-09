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
# Summary: ROSMAP meninges label transfer reference assembly
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

Stromal <- readRDS("stromal0129.rds")

DefaultAssay(Stromal) <- "RNA"

Stromal@meta.data$CellTypeBroad <- "Stromal"

Immune <- readRDS("MNGImmune0129.rds")

DefaultAssay(Immune) <- "RNA"

Immune@meta.data$CellTypeBroad <- "Immune"

Endos <- readRDS("endothelial0129.rds")

DefaultAssay(Endos) <- "RNA"

Endos@meta.data$CellTypeBroad <- "Endothelial"

metacols <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "cluster", "nCount_raw", "nFeature_raw",
              "id", "barcode", "batch", "projid", "subset", "class",
              "cell.type", "state", "sub.population", "ct", "CellTypeBroad")

s_merge <- list()

s_meta <-  list()

s_list <- c("Stromal" = Stromal, "Endos" = Endos, "Immune" = Immune)

names(s_list)

for(cell in names(s_list)){
  
  s_merge[[cell]] <- as.data.frame(GetAssayData(s_list[[cell]]))
  
  s_meta[[cell]] <- s_list[[cell]]@meta.data[,metacols]
  
}

s_merged <- list.cbind(s_merge)

s_metadata <- list.rbind(s_meta)

s <- CreateSeuratObject(counts = s_merged, meta.data = s_metadata)

s2 <- s %>% 
  NormalizeData() %>% 
  ScaleData() %>% 
  FindVariableFeatures()

s2 <- RunPCA(s2, features = s2@assays[["RNA"]]@meta.data[["var.features"]], npcs = 50)

s2 <- RunUMAP(s2, dims = 1:30)

SCP::CellDimPlot(s2, reduction = "umap", group.by = c("ct", "CellTypeBroad"),
                 theme_use = "theme_cowplot", legend.position = "right")



s3 <- s2

s3@assays$RNA <- as(object = s3[["RNA"]], Class = "Assay")

adata <- SCP::srt_to_adata(s3)

adata$write_h5ad("MNG_0129.h5ad")

saveRDS(s2, "MNG_0129.rds")