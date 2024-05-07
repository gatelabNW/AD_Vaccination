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
# Summary: Make spatial mixture plots to show Cell2Location results
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
  library(stardust)
  library(DescTools)
  library(ppcor)
  library(fastDummies)
  library(SPOTlight)
  library(gridExtra)
})

options(future.globals.maxSize=1048576000000)
load("helperfunctions.RData")
setwd("/path/to/C2L/output")

s <- readRDS("seurat_with_C2L_predictions.rds")

plts <- list()

for(donor in unique(s@meta.data$sample_id)){
  
  s2 <- subset(s, sample_id == donor)
  
  s2 <- subset(s2, manual_annotation != "")
  
  colnames(s2@meta.data)
  
  roughdeconv <- s2@meta.data[,39:56]
  
  roughdeconv$None <- ifelse(rowSums(roughdeconv) == 0, 1, 0)
  
  cnames <- colnames(roughdeconv)
  
  c2 <- as.matrix(GetTissueCoordinates(s2, image = paste0(donor)))
  
  colnames(c2) <- c("coord_y", "coord_x")
  
  if(donor == "AN1792.102.1"){ # this sample has negative coordinates?
    c2[,1] <- c2[,1] - min(c2[,1])
  }else{
    c2[,1] <- c2[,1]
  }
  
  c3 <- as.data.frame(c2)
  
  class(c2) <- "numeric"
  
  c3$layer <- s2@meta.data$manual_annotation
  
  colrs <- c("chartreuse", # Astros
             "darkgreen", # Endos
             "purple4", # FB
             "yellow", # interneurons
             "darkblue", "blue3", "blue", "dodgerblue2", "steelblue", "skyblue", "turquoise", "cyan", #ENs
             "red", #glia,
             "orange", "darkorange4",# OPC/Oligo
             "blueviolet", "mediumorchid", "magenta", # Immune
             "gray80" #None
  )
  
  fillcols <- c("gray20", "gray40", "gray60", "gray80", "gray90", "black", "white")
  
  plts[[donor]] <- plotSpatialScatterpie(
    x = c2,
    y = roughdeconv,
    cell_types = cnames,
    img = FALSE,
    pie_scale = 0.33) +
    geom_point(data = c3, mapping = aes(x = coord_x, y = coord_y, color = layer), size = 2, alpha = 0.25) +
    theme_cowplot() +
    theme(
      plot.background = element_rect(fill = "white"), 
      axis.line.x.bottom=element_line(color="black"),
      axis.line.y.left=element_line(color="black"),
      legend.text = element_text(color="black"),
      title = element_text(color="black"),
      axis.ticks = element_line(color="black"),
      axis.text = element_text(color="black")
    ) +
    scale_y_reverse() + 
    scale_fill_manual(values = colrs) + 
    scale_color_manual(values = fillcols) +
    labs(fill = "Cell Type", 
         x = "", 
         y = "") +
    ggtitle(paste0(donor, " Cell2Location Mixture Plot"))
  
  rm(s2)
  
  gc()
  
}

names(plts) <- gsub(".", "-", names(plts), fixed = TRUE)

pdf(file = paste0("Scatterpies.pdf"), width = 48, height = 48)

grid.arrange(plts$A34285, plts$A34717, plts$A34291A, plts$A34992A, plts$`A18-148`,plts$`AX21-92`, # ALL NNC
             plts$N35127N, plts$`A34933-2`, plts$A34995, plts$`A34583-2`, plts$A34644, plts$A35038, plts$`AN1792-102-15`, #nAD and Placebo
             plts$`AN1792-102-4`, plts$`AN1792-102-10`, plts$`AN1792-102-11`,  plts$`AN1792-102-17`,plts$`AN1792-102-19`, plts$`AN1792-102-22`, #lim
             plts$`AN1792-102-1`, plts$`AN1792-102-6`, plts$`AN1792-102-7`, plts$`AN1792-102-8`, plts$`AN1792-102-16`, plts$`AN1792-102-20`,
             plts$`AN1792-102-21`, nrow = 6, ncol = 5)

dev.off()