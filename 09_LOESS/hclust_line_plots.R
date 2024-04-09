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
# Summary: LOESS hierarchical clustering and cluster line plots (Ext Fig 2D-F)
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
# Load and scale LOESS predictions for each group

# Load and scale all predictions 
lim_merged <- read.csv(paste0(output_folder, "lim_predictions.csv"))
lim_scaled <- lim_merged
lim_scaled[,colnames(lim_scaled) != "amyloid"] <- scale(lim_scaled[,colnames(lim_scaled) != "amyloid"])
lim_scaled <- lim_scaled[,colSums(is.na(lim_scaled)) == 0]

ext_merged <- read.csv(paste0(output_folder, "ext_predictions.csv"))
ext_scaled <- ext_merged
ext_scaled[,colnames(ext_scaled) != "amyloid"] <- scale(ext_scaled[,colnames(ext_scaled) != "amyloid"])
ext_scaled <- ext_scaled[,colSums(is.na(ext_scaled)) == 0]

nAD_merged <- read.csv(paste0(output_folder, "nAD_predictions.csv"))
nAD_scaled <- nAD_merged
nAD_scaled[,colnames(nAD_scaled) != "amyloid"] <- scale(nAD_scaled[,colnames(nAD_scaled) != "amyloid"])
nAD_scaled <- nAD_scaled[,colSums(is.na(nAD_scaled)) == 0]

#-------------------------------------------------------------------------------
# Run hierarchical clustering and generate k = 2:10 clusters per group

# iAD-Limited
lim_dist <- dist(t(lim_scaled[,2:ncol(lim_scaled)]))
lim_clust <- hclust(lim_dist, method = "complete")

# Generate range of clusters
cluster_nums <- 2:10
lim_clust_list <- list()
for (k in cluster_nums) {
  cur_cut <- data.frame(cutree(lim_clust, k = k))
  lim_clust_list[[k]] <- cur_cut
}

lim_clust_merged <- as.data.frame(do.call(cbind, lim_clust_list[cluster_nums]))
colnames(lim_clust_merged) <- paste0("nclust_", cluster_nums)

# Save clusters
write.csv(lim_clust_merged, paste0(output_folder, "lim_clusters.csv"))

# iAD-Extensive
ext_dist <- dist(t(ext_scaled[,2:ncol(ext_scaled)]))
ext_clust <- hclust(ext_dist, method = "complete")

# Generate range of clusters
cluster_nums <- 2:10
ext_clust_list <- list()
for (k in cluster_nums) {
  cur_cut <- data.frame(cutree(ext_clust, k = k))
  ext_clust_list[[k]] <- cur_cut
}

ext_clust_merged <- as.data.frame(do.call(cbind, ext_clust_list[cluster_nums]))
colnames(ext_clust_merged) <- paste0("nclust_", cluster_nums)

# Save clusters
write.csv(ext_clust_merged, paste0(output_folder, "ext_clusters.csv"))

# nAD
nAD_dist <- dist(t(nAD_scaled[,2:ncol(nAD_scaled)]))
nAD_clust <- hclust(nAD_dist, method = "complete")

# Generate range of clusters
cluster_nums <- 2:10
nAD_clust_list <- list()
for (k in cluster_nums) {
  cur_cut <- data.frame(cutree(nAD_clust, k = k))
  nAD_clust_list[[k]] <- cur_cut
}

nAD_clust_merged <- as.data.frame(do.call(cbind, nAD_clust_list[cluster_nums]))
colnames(nAD_clust_merged) <- paste0("nclust_", cluster_nums)

# Save clusters
write.csv(nAD_clust_merged, paste0(output_folder, "nAD_clusters.csv"))

#-------------------------------------------------------------------------------
# Generate line plots for 6 clusters per group, using scaled predictions

# Define number of clusters to plot
n_clust <- 6

# Define cluster colors
colors <- hue_pal()(n_clust)

# Define nrow and ncol for ggarrange
nrow <- 2
ncol <- 3

# iAD-Limited
lim_clust_merged <- read.csv(paste0(output_folder, "lim_clusters.csv"), row.names = 1)
lim_long <- pivot_longer(lim_scaled, colnames(lim_scaled)[2:ncol(lim_scaled)], names_to = "gene", values_to = "expression")

# Define uniform y axis limits
min <- min(lim_long$expression)
max <- max(lim_long$expression)

# Generate line plots for each cluster, and plot average cluster expression
plt_list <- list()
for (clust in 1:n_clust) {
  cur_data <- lim_long[lim_long$gene %in% row.names(lim_clust_merged)[lim_clust_merged[,paste0("nclust_", n_clust)] == clust],]
  
  if (length(unique(cur_data$gene)) > 1) {
    avg_data <- cur_data %>% group_by(amyloid) %>% summarize(average = mean(expression))
    
    plt <- ggplot(cur_data, aes(x = amyloid, y = expression)) +
      theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
      geom_line(linewidth = 0.3, aes(group = gene), color = colors[clust]) +
      geom_line(data = avg_data,
                aes(x = amyloid, y = average),
                stat="smooth", method = "loess", span = 0.75, se = FALSE,
                color = darken(colors[[clust]], amount = 0.3), linewidth = 2) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
            axis.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
            axis.text.y = element_text(size = 15)) +
      ggtitle(paste0("Cluster ", clust, ": ", length(unique(cur_data$gene)), " Genes")) + theme(aspect.ratio = 1) + ylim(min, max)
  } else {
    plt <- ggplot(cur_data, aes(x = amyloid, y = expression)) +
      theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
      geom_line(linewidth = 0.3, aes(group = gene), color = colors[clust]) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
            axis.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
            axis.text.y = element_text(size = 15)) +
      ggtitle(paste0("Cluster ", clust, ": ", length(unique(cur_data$gene)), " Genes")) + theme(aspect.ratio = 1) + ylim(min, max)
  }
  plt_list[[clust]] <- plt
}

# Ext Fig 2E
pdf(paste0(output_folder, "lim_plots/lim_cluster_lineplots_", n_clust, "_same_scale.pdf"), height = 10*nrow, width =10*ncol)
ggarrange(plotlist = plt_list, ncol = ncol, nrow = nrow)
dev.off()

# iAD-Extensive
ext_clust_merged <- read.csv(paste0(output_folder, "ext_clusters.csv"), row.names = 1)
ext_long <- pivot_longer(ext_scaled, colnames(ext_scaled)[2:ncol(ext_scaled)], names_to = "gene", values_to = "expression")

# Define uniform y axis limits
min <- min(ext_long$expression)
max <- max(ext_long$expression)

# Generate line plots for each cluster, and plot average cluster expression
plt_list <- list()
for (clust in 1:n_clust) {
  cur_data <- ext_long[ext_long$gene %in% row.names(ext_clust_merged)[ext_clust_merged[,paste0("nclust_", n_clust)] == clust],]
  
  if (length(unique(cur_data$gene)) > 1) {
    avg_data <- cur_data %>% group_by(amyloid) %>% summarize(average = mean(expression))
    
    plt <- ggplot(cur_data, aes(x = amyloid, y = expression)) +
      theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
      geom_line(linewidth = 0.3, aes(group = gene), color = colors[clust]) +
      geom_line(data = avg_data,
                aes(x = amyloid, y = average),
                stat="smooth", method = "loess", span = 0.75, se = FALSE,
                color = darken(colors[[clust]], amount = 0.3), linewidth = 2) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
            axis.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
            axis.text.y = element_text(size = 15)) +
      ggtitle(paste0("Cluster ", clust, ": ", length(unique(cur_data$gene)), " Genes")) + theme(aspect.ratio = 1) + ylim(min, max)
  } else {
    plt <- ggplot(cur_data, aes(x = amyloid, y = expression)) +
      theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
      geom_line(linewidth = 0.3, aes(group = gene), color = colors[clust]) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
            axis.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
            axis.text.y = element_text(size = 15)) +
      ggtitle(paste0("Cluster ", clust, ": ", length(unique(cur_data$gene)), " Genes")) + theme(aspect.ratio = 1) + ylim(min, max)
  }
  plt_list[[clust]] <- plt
}

# Ext Fig 2F
pdf(paste0(output_folder, "ext_plots/ext_cluster_lineplots_", n_clust, "_same_scale.pdf"), height = 10*nrow, width =10*ncol)
ggarrange(plotlist = plt_list, ncol = ncol, nrow = nrow)
dev.off()

# nAD
nAD_clust_merged <- read.csv(paste0(output_folder, "nAD_clusters.csv"), row.names = 1)
nAD_long <- pivot_longer(nAD_scaled, colnames(nAD_scaled)[2:ncol(nAD_scaled)], names_to = "gene", values_to = "expression")

# Define uniform y axis limits
min <- min(nAD_long$expression)
max <- max(nAD_long$expression)

# Generate line plots for each cluster, and plot average cluster expression
plt_list <- list()
for (clust in 1:n_clust) {
  cur_data <- nAD_long[nAD_long$gene %in% row.names(nAD_clust_merged)[nAD_clust_merged[,paste0("nclust_", n_clust)] == clust],]
  
  if (length(unique(cur_data$gene)) > 1) {
    avg_data <- cur_data %>% group_by(amyloid) %>% summarize(average = mean(expression))
    
    plt <- ggplot(cur_data, aes(x = amyloid, y = expression)) +
      theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
      geom_line(linewidth = 0.3, aes(group = gene), color = colors[clust]) +
      geom_line(data = avg_data,
                aes(x = amyloid, y = average),
                stat="smooth", method = "loess", span = 0.75, se = FALSE,
                color = darken(colors[[clust]], amount = 0.3), linewidth = 2) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
            axis.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
            axis.text.y = element_text(size = 15)) +
      ggtitle(paste0("Cluster ", clust, ": ", length(unique(cur_data$gene)), " Genes")) + theme(aspect.ratio = 1) + ylim(min, max)
  } else {
    plt <- ggplot(cur_data, aes(x = amyloid, y = expression)) +
      theme(panel.background = element_rect(fill = 'white'), axis.line = element_line(color = 'black')) +
      geom_line(linewidth = 0.3, aes(group = gene), color = colors[clust]) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
            axis.text.x = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
            axis.text.y = element_text(size = 15)) +
      ggtitle(paste0("Cluster ", clust, ": ", length(unique(cur_data$gene)), " Genes")) + theme(aspect.ratio = 1) + ylim(min, max)
  }
  plt_list[[clust]] <- plt
}

# Ext Fig 2D
pdf(paste0(output_folder, "nAD_plots/nAD_cluster_lineplots_", n_clust, "_same_scale.pdf"), height = 10*nrow, width =10*ncol)
ggarrange(plotlist = plt_list, ncol = ncol, nrow = nrow)
dev.off()









