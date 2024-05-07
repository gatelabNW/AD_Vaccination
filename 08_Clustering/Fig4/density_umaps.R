# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                         AN1792 Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 04-04-2024
# Written by: Anne Forsyth
# Summary: Density UMAPs for Lecanemab and CAA control microglia (Fig 3K)
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages ({
  library('Seurat')
  library('glue')
  library('dplyr')
  library('stringr')
  library("Matrix")
  library("data.table")
  library('sctransform')
  library('ggplot2')
  library("randomcoloR")
  library("ggpubr")
  library("tidyverse")
  library("introdataviz")
  library("scales")
  library("ggplotify")
})

output_dir <- "/path/to/general/output/folder/"

# Load Seurat object
load("/path/to/s_microglia_new_clean")

# Define function to calculate density (https://slowkow.com/notes/ggplot2-color-by-density/)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Define LCMB vs. CAA variable
s$sample_short <- str_split_fixed(s$sample_id, "\\.", 3)[,3]
s@meta.data <- s@meta.data %>% mutate(group = case_when(sample_short %in% c("A1", "A3", "A4", "A9") ~ "CAA",
                                                        sample_short %in% c("B1", "B3", "B4", "B9") ~ "LCMB"))


# Lecanemab UMAP 
temp_s <- subset(s, group == "LCMB")
data <- temp_s@reductions$integratedUMAP@cell.embeddings %>% as.data.frame
data$density <- get_density(data$integratedUMAP_1, data$integratedUMAP_2, n = 100)
plt <- ggplot(data, aes(x = integratedUMAP_1, y = integratedUMAP_2)) + 
  geom_point(aes(fill = density), alpha = 0.5, shape = 21, stroke = NA, size = 4) +
  stat_density_2d(geom = "density_2d", linewidth = 0.5, color = "gray30") +
  theme_classic() + viridis::scale_fill_viridis(option = "magma") +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(aspect.ratio = 1)

pdf(paste0(output_dir, "microglia_clusters_new_clean/LCMB_density_umap.pdf"), height = 10, width = 10)
print(plt)
dev.off()

# CAA control UMAP
temp_s <- subset(s, group == "CAA")
data <- temp_s@reductions$integratedUMAP@cell.embeddings %>% as.data.frame
data$density <- get_density(data$integratedUMAP_1, data$integratedUMAP_2, n = 100)
plt <- ggplot(data, aes(x = integratedUMAP_1, y = integratedUMAP_2)) + 
  geom_point(aes(fill = density), alpha = 0.5, shape = 21, stroke = NA, size = 4) +
  stat_density_2d(geom = "density_2d", linewidth = 0.5, color = "gray30") +
  theme_classic() + viridis::scale_fill_viridis(option = "magma") +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(aspect.ratio = 1)

pdf(paste0(output_dir, "microglia_clusters_new_clean/CAA_density_umap.pdf"), height = 10, width = 10)
print(plt)
dev.off()