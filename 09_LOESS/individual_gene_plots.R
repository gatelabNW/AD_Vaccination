# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                         AN1792 Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 04-03-2024
# Written by: Anne Forsyth
# Summary: LOESS line plots with confidence interval for gene expression over amyloid density (Fig 2K)
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("ggplot2")
  library("gridGraphics")
  library("ggplotify")
  library("ggpubr")
  library('ggnewscale')
  library("plotrix")
  library("stringi")
})

# Load integrated Seurat object 
s <- readRDS("/path/to/all_samples_03.rds")

# Define filter operator
`%notin%` <- Negate(`%in%`)

# Define output folder for amyloid density plots
output_folder <- "/path/to/graymatter/amyloiddensity/loessplots/"

#-------------------------------------------------------------------------------
# Generate LOESS line plots for genes of interest, using amyloid density

# Remove amyloid "exclude" spots and subset for gray matter 
s <- subset(s, amyloid_filter == "include" & manual_annotation %notin% c("white", "meninges"))

# Recorrect SCT counts across samples (layers are split by sample ID)
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}
s <- PrepSCTFindMarkers(object = s)

# Extract SCT data
all_data <- GetAssayData(s, assay = "SCT", layer = "data")

# Define genes to plot
plt_genes <- c("AQP4", "C3", "CD74", "SPARC", "APOE", "TREM2", "CD68", "BTK", "CGAS", "TGFBR1", "CXCL8", "PYCARD", "PROS1", "ITGAX", "TYROBP", "C3AR1", "IFNGR1", "MERTK")

# Format data for LOESS
plt_data <- all_data[plt_genes, ]
plt_data <- data.frame(plt_data)
plt_data <- t(plt_data)
plt_data <- data.frame(plt_data)
plt_data$row_name <- stri_replace_last_fixed(row.names(plt_data), ".", "-")
row.names(plt_data) <- plt_data$row_name
sum(row.names(plt_data) != row.names(s@meta.data))
sum(is.na(plt_data))
plt_data$amyloid <- s@meta.data$amyloid_fluo
plt_data$condition <- s@meta.data$condition_clearance

# Exclude NNC and placebo
plt_data <- plt_data[plt_data$condition %notin% c("NNC", "Placebo"),]

# Format condition variable for plot names
plt_data <- plt_data %>% mutate(condition = case_when(condition == "nAD" ~ "nAD",
                                                      condition == "lim" ~ "iAD-Limited",
                                                      condition == "ext" ~ "iAD-Extensive"))
plt_data$condition <- factor(plt_data$condition, levels = c("nAD", "iAD-Limited", "iAD-Extensive"))

# Define plot colors
colors <- c("royalblue4", "red4", "darkorange3")
names(colors) <- levels(plt_data$condition)

# Generate LOESS line plots with confidence interval (Fig 2K)
for (gene in plt_genes) {
  plt_data$gene <- plt_data[,gene]
  
  plt <- ggplot(data = plt_data, mapping = aes(x = amyloid, y = gene, group = condition, color = condition)) +
    stat_smooth(method = "loess", se = TRUE, mapping = aes(fill = condition)) +
    scale_fill_manual(values = colors) + scale_color_manual(values = colors) + 
    labs(y = "SCTransform Expression", x = "Amyloid Density", color = NULL) +
    theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black'),
          plot.title = element_text(hjust = 0.5, size = 15),
          axis.text = element_text(size = 12), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
          legend.text = element_text(size = 12), aspect.ratio = 1) + ggtitle(gene) + guides(fill = "none") +
    geom_vline(xintercept = 153, linetype = 3)
  
  pdf(paste0(output_folder, gene, ".pdf"), height = 10, width = 10)
  print(plt)
  dev.off()
}

