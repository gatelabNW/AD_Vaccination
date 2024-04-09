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
# Summary: UpSet plot and heatmap of positive LOESS clusters (Fig 2J, Ext Fig 2G)
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggrepel")
  library("ggthemes")
  library("ggpubr")
  library("Seurat")
  library("factoextra")
  library("scales")
  library("patchwork")
  library("ggnewscale")
  library("clusterProfiler")
  library("pheatmap")
  library("colorspace")
  library("stringi")
  library("randomcoloR")
  library("UpSetR")
})

# Define output folder
output_folder <- "/path/to/loess/output/folder/"

# Define filter operator
`%notin%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
# Load cluster data and generate UpSet plot for clusters with positive trend

# Load nAD clusters, filter for clusters 1, 3
nAD_clusters <- read.csv(paste0(output_folder, "nAD_clusters.csv"), row.names = 1)
nAD_clusters$gene <- row.names(nAD_clusters)
nAD_clusters <- nAD_clusters[nAD_clusters$nclust_6 %in% c(1, 3),]

# Load iAD-Limited clusters, filter for clusters 1, 3, 4
lim_clusters <- read.csv(paste0(output_folder, "lim_clusters.csv"), row.names = 1)
lim_clusters$gene <- row.names(lim_clusters)
lim_clusters <- lim_clusters[lim_clusters$nclust_6 %in% c(1, 3, 4),]

# Load iAD-Extensive clusters, filter for cluster 1
ext_clusters <- read.csv(paste0(output_folder, "ext_clusters.csv"), row.names = 1)
ext_clusters$gene <- row.names(ext_clusters)
ext_clusters <- ext_clusters[ext_clusters$nclust_6 == 1,]

# Combine cluster data in list for UpSet plot
clust_data <- list(nAD = nAD_clusters, lim = lim_clusters, ext = ext_clusters)
plt_sets <- list()
for (group in names(clust_data)) {
  cur_data <- clust_data[[group]]
  for (clust in unique(cur_data$nclust_6)) {
    plt_sets[[paste0(group, "_", clust)]] <- unique(cur_data$gene[cur_data$nclust_6 == clust])
  }
}

# Define bar colors
set.seed(1)
bar_colors <- distinctColorPalette(23)

# UpSet plot of all genes per cluster 
plt <- UpSetR::upset(fromList(plt_sets), order.by = "freq", nsets = length(plt_sets), nintersects = 23, 
                     sets.bar.color = "darkblue", main.bar.color = bar_colors)

# Ext Fig 2G
pdf(paste0(output_folder, "positive_cluster_comparison/upset_6_clusters.pdf"), height = 10, width = 15)
print(plt)
dev.off()

#-------------------------------------------------------------------------------
# Load and scale predictions for each group

# iAD-Limited
lim_merged <- read.csv(paste0(output_folder, "lim_predictions.csv"))
lim_scaled <- lim_merged
lim_scaled[,colnames(lim_scaled) != "amyloid"] <- scale(lim_scaled[,colnames(lim_scaled) != "amyloid"])
lim_scaled <- lim_scaled[,colSums(is.na(lim_scaled)) == 0]

# iAD-Extensive
ext_merged <- read.csv(paste0(output_folder, "ext_predictions.csv"))
ext_scaled <- ext_merged
ext_scaled[,colnames(ext_scaled) != "amyloid"] <- scale(ext_scaled[,colnames(ext_scaled) != "amyloid"])
ext_scaled <- ext_scaled[,colSums(is.na(ext_scaled)) == 0]

# nAD
nAD_merged <- read.csv(paste0(output_folder, "nAD_predictions.csv"))
nAD_scaled <- nAD_merged
nAD_scaled[,colnames(nAD_scaled) != "amyloid"] <- scale(nAD_scaled[,colnames(nAD_scaled) != "amyloid"])
nAD_scaled <- nAD_scaled[,colSums(is.na(nAD_scaled)) == 0]

# Combine gene expression columns 
lim_long <- pivot_longer(lim_scaled, colnames(lim_scaled)[2:ncol(lim_scaled)], names_to = "gene", values_to = "expression")
ext_long <- pivot_longer(ext_scaled, colnames(ext_scaled)[2:ncol(ext_scaled)], names_to = "gene", values_to = "expression")
nAD_long <- pivot_longer(nAD_scaled, colnames(nAD_scaled)[2:ncol(nAD_scaled)], names_to = "gene", values_to = "expression")

# Find unique genes per cluster
unique_genes <- list()
for (cluster in names(plt_sets)) {
  cur_genes <- plt_sets[[cluster]]
  other_genes <- unlist(plt_sets[names(plt_sets) != cluster])
  unique_genes[[cluster]] <- cur_genes[cur_genes %notin% other_genes]
  print(cluster)
  print(length(unique_genes[[cluster]]))
}

# Exclude genes with NA scaled expression in any group from plots
shared_genes <- intersect(unique(lim_long$gene), unique(ext_long$gene))
shared_genes <- intersect(shared_genes, unique(nAD_long$gene))
ordered <- unlist(unique_genes)
ordered <- ordered[ordered %in% shared_genes]

# Filter for genes in unique sets
lim_long <- lim_long[lim_long$gene %in% ordered,]
ext_long <- ext_long[ext_long$gene %in% ordered,]
nAD_long <- nAD_long[nAD_long$gene %in% ordered,]

# Order genes by unique cluster sets
lim_long$gene <- factor(lim_long$gene, levels = ordered)
lim_long <- lim_long %>% arrange(gene)
nAD_long$gene <- factor(nAD_long$gene, levels = ordered)
nAD_long <- nAD_long %>% arrange(gene)
ext_long$gene <- factor(ext_long$gene, levels = ordered)
ext_long <- ext_long %>% arrange(gene)

# Filter amyloid density to minimum total range (for uniform x axis limits)
amyloid_thresh <- min(c(max(lim_long$amyloid), max(ext_long$amyloid), max(nAD_long$amyloid)))
lim_long <- lim_long[lim_long$amyloid <= amyloid_thresh,]
ext_long <- ext_long[ext_long$amyloid <= amyloid_thresh,]
nAD_long <- nAD_long[nAD_long$amyloid <= amyloid_thresh,]

# Define limits for uniform gene expression scale
min <- min(c(nAD_long$expression, lim_long$expression, ext_long$expression))
max <- max(c(nAD_long$expression, lim_long$expression, ext_long$expression))

# Add index variable for numeric cutoffs
lim_long$index <- 1:nrow(lim_long)
ext_long$index <- 1:nrow(ext_long)
nAD_long$index <- 1:nrow(nAD_long)

# Define y axis cutoffs for unique gene sets (this will be the same for all groups)
cutoffs <- c()
for (clust in names(unique_genes)) {
  cutoffs <- c(cutoffs, lim_long$gene[max(lim_long$index[lim_long$gene %in% unique_genes[[clust]]])])
}

# Define y axis labels for unique gene sets
ticks <- list()
for (i in 1:length(cutoffs)) {
  if (i == 1) {
    ticks[[names(unique_genes)[i]]] <- cutoffs[i]/2
  } else {
    ticks[[names(unique_genes)[i]]] <- cutoffs[i] - (cutoffs[i] - cutoffs[i - 1])/2
  }
}

# Generate heatmaps for each group
lim_heatmap <- ggplot(lim_long, aes(x = amyloid, y = gene, fill = expression)) +
  theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
  geom_raster(interpolate=TRUE) +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 15, vjust = -1), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12)) +
  scale_fill_gradientn(limits = c(min, max), colours=c("turquoise4", "paleturquoise", "lightyellow", "plum1", "mediumorchid3")) +
  ggtitle(paste0("iAD-Limited: " , length(unique(lim_long$gene)), " Genes")) + 
  labs(x = "Amyloid Density", y ="Gene", fill = "Predicted Expression") + 
  theme(aspect.ratio = 1) + geom_hline(yintercept = cutoffs, linetype = 3) + scale_x_continuous(expand = c(0,0))

ext_heatmap <- ggplot(ext_long, aes(x = amyloid, y = gene, fill = expression)) +
  theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
  geom_raster(interpolate=TRUE) +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 15, vjust = -1), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12)) +
  scale_fill_gradientn(limits = c(min, max), colours=c("turquoise4", "paleturquoise", "lightyellow", "plum1", "mediumorchid3")) +
  ggtitle(paste0("iAD-Extensive: " , length(unique(ext_long$gene)), " Genes")) + 
  labs(x = "Amyloid Density", y ="Gene", fill = "Predicted Expression") + 
  theme(aspect.ratio = 1) + geom_hline(yintercept = cutoffs, linetype = 3) + scale_x_continuous(expand = c(0,0))

nAD_heatmap <- ggplot(nAD_long, aes(x = amyloid, y = gene, fill = expression)) +
  theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
  geom_raster(interpolate=TRUE) +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 15, vjust = -1), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12)) +
  scale_fill_gradientn(limits = c(min, max), colours=c("turquoise4", "paleturquoise", "lightyellow", "plum1", "mediumorchid3")) +
  ggtitle(paste0("nAD: " , length(unique(nAD_long$gene)), " Genes")) + 
  labs(x = "Amyloid Density", y ="Gene", fill = "Predicted Expression") + 
  theme(aspect.ratio = 1) + geom_hline(yintercept = cutoffs, linetype = 3) +
  annotate("text", x = -10, y = unlist(ticks), label = names(ticks)) + scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(0, max(nAD_long$amyloid)),  clip = 'off') 

# Fig 2J
pdf(paste0(output_folder, "positive_cluster_comparison/heatmaps_edited.pdf"), height = 10, width = 35)
print((nAD_heatmap + theme(legend.position = "none", axis.title.x = element_blank())) + 
        (lim_heatmap + theme(legend.position = "none", axis.text.y = element_blank())) + 
        (ext_heatmap + theme(axis.title.x = element_blank(), axis.text.y = element_blank())))
dev.off()

