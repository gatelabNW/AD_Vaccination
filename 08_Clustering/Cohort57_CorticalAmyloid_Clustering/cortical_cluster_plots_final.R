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
# Summary: Cortical WNN cluster marker plots (Fig 4K, Fig 4H, Ext Fig 4J, Ext Fig 4N)
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
  library("randomcoloR")
  library("ggpubr")
  library("gridExtra")
  library("ggplot2")
  library("ggnewscale")
  library("ggheatmap")
  library("radiant.data")
  library("ggtree")
  library("aplot")
  library("ggrepel")
  library("ggplotify")
  library("grid")
  library("ggdendroplot")
})

# Define input and output folders
marker_dir <- "/path/to/findmarkers/results/"
output_dir <- "/path/to/cortical/wnn/plots/"

# Load Seurat object with WNN
s <- readRDS("/path/to/plaque_wnn.rds")

# Format cluster names
s$cluster <- paste0("plaque_", s$wsnn_res.0.5)
s$cluster <- factor(s$cluster, levels = paste0("plaque_", 0:14))
Idents(s) <- "cluster"

# For PrepSCTFindMarkers
for (name in names(s@assays$SCT@SCTModel.list)) {
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}

# Load significant marker results (calculated using all clusters)
all_sig_markers <- read.csv(paste0(marker_dir, "all_sig_markers.csv"))

#-------------------------------------------------------------------------------
# Violin plot of modality weights

# Extract modality weight meta data
data <- s@meta.data[,c("SCT.weight", "Protein.weight", "cluster")]
weight <- c(data$SCT.weight, data$Protein.weight)
cluster <- rep(data$cluster, 2)
modality <- c(rep("RNA", nrow(data)), rep("Protein", nrow(data)))
plt_data <- data.frame(Weight = weight, Cluster = cluster, Modality = modality)
plt <- ggplot(plt_data, aes(x = Cluster, y = Weight, fill = Modality)) + geom_split_violin() + scale_fill_manual(values = c("lightseagreen", "orchid3")) + 
  theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black'), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Plaque WNN Cluster Modality Weights: Resolution = 0.5") + theme(plot.title = element_text(hjust = 0.5))

# Ext Fig 4I
pdf(paste0(output_dir, "plaque_cluster_modality_weights.pdf"), height = 10, width = 15)
print(plt)
dev.off()

#-------------------------------------------------------------------------------
# Density UMAP for Lecanemab and CAA control (Fig 4H)

# Define function to calculate density
# https://slowkow.com/notes/ggplot2-color-by-density/
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Define LCMB vs. CAA variable
s$sample_short <- str_split_fixed(s$sample_id, "\\.", 2)[,2]
s@meta.data <- s@meta.data %>% mutate(group = case_when(sample_short %in% c("A1", "A3", "A4", "A9") ~ "CAA",
                                                        sample_short %in% c("B1", "B3", "B4", "B9") ~ "LCMB"))

# Lecanemab
temp_s <- subset(s, group == "LCMB")
data <- temp_s@reductions$wnn.umap@cell.embeddings %>% as.data.frame
data$density <- get_density(data$wnnUMAP_1, data$wnnUMAP_2, n = 100)

plt <- ggplot(data, aes(x = wnnUMAP_1, y = wnnUMAP_2)) + 
  geom_point(aes(fill = density), alpha = 0.5, shape = 21, stroke = NA, size = 4) +
  stat_density_2d(geom = "density_2d", linewidth = 0.5, color = "gray30") +
  theme_classic() + viridis::scale_fill_viridis(option = "magma") +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(aspect.ratio = 1)

pdf(paste0(output_dir, "LCMB_density_umap.pdf"), height = 10, width = 10)
print(plt)
dev.off()

# CAA control 
temp_s <- subset(s, group == "CAA")
data <- temp_s@reductions$wnn.umap@cell.embeddings %>% as.data.frame
data$density <- get_density(data$wnnUMAP_1, data$wnnUMAP_2, n = 100)

plt <- ggplot(data, aes(x = wnnUMAP_1, y = wnnUMAP_2)) + 
  geom_point(aes(fill = density), alpha = 0.5, shape = 21, stroke = NA, size = 4) +
  stat_density_2d(geom = "density_2d", linewidth = 0.5, color = "gray30") +
  theme_classic() + viridis::scale_fill_viridis(option = "magma") +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(aspect.ratio = 1)

pdf(paste0(output_dir, "CAA_density_umap.pdf"), height = 10, width = 10)
print(plt)
dev.off()

#-------------------------------------------------------------------------------
# Generate heatmap for subset of clusters

# Subset for clusters of interest
temp <- subset(s, wsnn_res.0.5 %in% c(5, 6, 7, 10, 11, 12, 13))

# Recorrect SCT counts across samples
temp <- PrepSCTFindMarkers(temp)

# Identify top 3 marker genes per cluster based on PFC
t3_genes_to_plot <- c()
for(clust in unique(temp$wsnn_res.0.5)){
  print(clust)
  data <- dplyr::filter(all_sig_markers, cluster==clust)
  print(sum(data$PFC == Inf))
  data$PFC[data$PFC == Inf] <- (-log10(min(data$BH[data$BH != 0])) + 10)*abs(data$avg_log2FC[data$PFC == Inf])
  print(sum(data$PFC == Inf))
  data <- data %>% dplyr::arrange(desc(PFC))
  t3_genes_to_plot <- c(t3_genes_to_plot, data$gene|>head(3))
}
t3_genes_to_plot <- t3_genes_to_plot|>unique()

# Heatmap with dendrogram, using top 3 markers per cluster
data <- GetAssayData(temp, assay = "SCT", layer = "data")
data <- as.matrix(data)
sum(colnames(data) != row.names(temp@meta.data))
colnames(data) <- temp@meta.data$cluster
aggregated <- t(rowsum(t(data), group = colnames(data)))
for (col in colnames(aggregated)) {
  aggregated[,col] <- aggregated[,col] / sum(colnames(data) == col)
} 

# Subset for marker genes and scale gene expression across clusters
plt_data <- aggregated[row.names(aggregated) %in% t3_genes_to_plot,]
scaled <- t(plt_data)
scaled <- scale(scaled)
scaled <- t(scaled)

# Format data for ggplot
data <- as.data.frame(scaled) %>% rownames_to_column(var = "gene") %>% tidyr::gather(key = "cluster",value = "expression",-gene)

# Cluster rows
row_clust <- hclust(dist(as.matrix(scaled)))
roworder <- row_clust$labels[row_clust$order]
data$gene <- factor(data$gene, levels = roworder)

# Cluster columns
col_clust <- hclust(dist(as.matrix(t(scaled))))
colorder <- col_clust$labels[col_clust$order]
data$cluster <- factor(data$cluster, levels = colorder)

# Generate continuous x and y values (required to add dendrogram)
data$y <- NA
i <- 1
for (gene in levels(data$gene)) {
  data$y[data$gene == gene] <- i
  i <- i + 1
}

data$x <- NA
i <- 1
for (cluster in levels(data$cluster)) {
  data$x[data$cluster == cluster] <- i
  i <- i + 1
}

# Initialize plot
plt <- ggplot(data, mapping = aes(x = x, y = y)) + geom_tile(aes(fill = expression), color = "white", linewidth = 0.5) +
  scale_fill_gradientn(colours=c("#330099", "#99CCFF", "#FFFFCC", "#FFCC66", "#CC0000"))

# Add dendrograms
plt <- plt + geom_dendro(row_clust, pointing = "side", xlim=c(7.5, 9.5))
plt <- plt + geom_dendro(col_clust, pointing = "updown", ylim=c(21.5, 22))

# Edit plot
plt <- plt + labs(x = "Plaque WNN Cluster", y = "Marker Gene") +
  theme(panel.background = element_rect(fill = 'white')) + 
  ggtitle("Scaled Average SCT Expression of FindMarker Genes for WNN Clusters") + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 2), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_equal() + labs(fill = "Scaled Expression")

# Fig 4K
ggsave(filename = "plaque_dendro_heatmap.pdf", plot = plt, device = "pdf", path = output_dir, scale = 0.5, height = 20, width = 15)

# Corresponding IBA1 heatmap 
data <- temp@meta.data
data$cluster <- factor(data$cluster, levels = colorder)
summary <- data %>% dplyr::group_by(cluster) %>% dplyr::summarize(iba1 = mean(iba1_fluo))

plt <- ggplot(summary, aes(x = cluster, y = 1)) + geom_tile(aes(fill = iba1), color = "white", linewidth = 0.5) +
  scale_fill_gradientn(colours=c("#330099", "#99CCFF", "#FFFFCC", "#FFCC66", "#CC0000")) +
  theme(panel.background = element_rect(fill = 'white'), plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) +
  coord_equal() + labs(fill = "Average IBA1 Density")

# Ext Fig 4N
ggsave(filename = "plaque_dendro_heatmap_iba1.pdf", plot = plt, device = "pdf", path = output_dir, scale = 0.5, height = 5, width = 15)

#-------------------------------------------------------------------------------
# Generate dot plot for all clusters, top 2 marker genes by PFC 

# Recorrect SCT counts across samples
s <- PrepSCTFindMarkers(s)

# Get top 2 genes from every cluster by PFC 
t2_genes_to_plot <- c()
for(clust in all_sig_markers$cluster|>unique()){
  data <- dplyr::filter(all_sig_markers, cluster==clust)
  print(sum(data$PFC == Inf))
  data$PFC[data$PFC == Inf] <- (-log10(min(data$BH[data$BH != 0])) + 10)*abs(data$avg_log2FC[data$PFC == Inf])
  print(sum(data$PFC == Inf))
  data <- data %>% dplyr::arrange(desc(PFC))
  t2_genes_to_plot <- c(t2_genes_to_plot, data$gene|>head(2))
}
t2_genes_to_plot <- t2_genes_to_plot|>unique()

# Ext Fig 4J
dotplot_out <- paste0(output_dir, "dot_plot_pfc_genes_top_2.pdf")
dot_plt <- DotPlot(object = s, features = t2_genes_to_plot) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

pdf(dotplot_out, width = 10, height = 8)
plot(dot_plt)
dev.off()



