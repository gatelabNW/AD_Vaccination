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
# Summary: Merge microglia clusters (Fig 3M, Ext Fig 3J)
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
  library("ggdendroplot")
})

# Define output folder
output_dir <- "/path/to/general/output/folder/"

# Load clean microglia cluster object
load("/path/to/s_microglia_new_clean")

#-------------------------------------------------------------------------------
# Merge clusters 

# Merge clusters 1, 5, 8, 11
s$merged_cluster <- as.character(s$microglia_cluster)
s$merged_cluster[s$microglia_cluster %in% paste0("microglia_", c(1, 5, 8, 11))] <- "microglia_merged"

# Save Seurat object with updated metadata
save(s, file = "/path/to/s_microglia_new_clean")

#-------------------------------------------------------------------------------
# Find marker genes for merged clusters

# Define output folder, recorrect SCT counts across layers, set idents
marker_dir <- paste0(output_dir, "microglia_clusters_new_clean/merged_markers/")
s <- PrepSCTFindMarkers(s)
Idents(s) <- "merged_cluster"

# Run FindMarkers
markers <- list()
sig_markers <- list()
for (clust in unique(s@meta.data$merged_cluster)){
  print(clust)
  results <- FindMarkers(s, ident.1 = clust, min.pct = 0.25, recorrect_umi = FALSE, only.pos = TRUE)
  results$BH <- p.adjust(results$p_val, method = "BH")
  results$cluster <- clust
  results$gene <- rownames(results)
  results$PFC <- -log10(results$BH)*abs(results$avg_log2FC)
  sig_results <- dplyr::filter(results, BH < 0.05)
  markers[[clust]] <- results
  sig_markers[[clust]] <- sig_results
}

# Compile results
all_sig_markers <- data.table::rbindlist(sig_markers)|>as.data.frame()
all_markers <- data.table::rbindlist(markers)|>as.data.frame()

# Define output file names
sig_out_file <- paste0(marker_dir, "all_sig_markers.csv")
out_file <- paste0(marker_dir, "all_markers.csv")

# Remove contamination
all_sig_markers_filtered <- all_sig_markers[!grepl("^MT-", all_sig_markers$gene),]
all_sig_markers_filtered <- all_sig_markers_filtered[!grepl("^HB", all_sig_markers_filtered$gene),]

# Save results
write.csv(all_sig_markers_filtered, sig_out_file, row.names = FALSE, quote = FALSE)
write.csv(all_markers, out_file, row.names = FALSE, quote = FALSE)

# Load significant marker results
all_sig_markers <- read.csv(paste0(marker_dir, "all_sig_markers.csv"))

# Identify top 5 markers per cluster based on PFC
marker_genes_to_plot <- c()
for(clust in all_sig_markers$cluster|>unique()){
  data <- dplyr::filter(all_sig_markers, cluster==clust)
  
  print(sum(data$PFC == Inf))
  data$PFC[data$PFC == Inf] <- (-log10(min(data$BH[data$BH != 0])) + 10)*abs(data$avg_log2FC[data$PFC == Inf])
  print(sum(data$PFC == Inf))
  data <- data %>% dplyr::arrange(desc(PFC))
  
  marker_genes_to_plot <- c(marker_genes_to_plot, data$gene|>head(5))
}

# Also plot homeostatic markers
marker_genes_to_plot <- c(marker_genes_to_plot, c("CSF1R", "P2RY12", "CX3CR1"))
marker_genes_to_plot <- marker_genes_to_plot|>unique()

# Seurat dot plot (Ext Fig 3J)
dotplot_out <- paste0(marker_dir, "dot_plot.pdf")
dot_plt <- DotPlot(object = s, features = marker_genes_to_plot) + theme(axis.text.x = element_text(angle = 90, size = 8))

pdf(dotplot_out, width = 15, height = 8)
plot(dot_plt)
dev.off()

#-------------------------------------------------------------------------------
# Heatmap for specific clusters

# Subset for clusters 2, 3, 4, 6, recorrect SCT counts across layers 
s <- subset(s, merged_cluster %in% paste0("microglia_", c(2:4, 6)))
s <- PrepSCTFindMarkers(s)

# Get top 5 markers per cluster based on PFC
marker_genes_to_plot <- c()
for(clust in paste0("microglia_", c(2:4, 6))){
  data <- dplyr::filter(all_sig_markers, cluster==clust)
  
  print(sum(data$PFC == Inf))
  data$PFC[data$PFC == Inf] <- (-log10(min(data$BH[data$BH != 0])) + 10)*abs(data$avg_log2FC[data$PFC == Inf])
  print(sum(data$PFC == Inf))
  data <- data %>% dplyr::arrange(desc(PFC))
  
  marker_genes_to_plot <- c(marker_genes_to_plot, data$gene|>head(5))
}

marker_genes_to_plot <- marker_genes_to_plot|>unique()

# Aggregated heatmap using average gene expression
data <- GetAssayData(s, assay = "SCT", layer = "data")
data <- as.matrix(data)
sum(colnames(data) != row.names(s@meta.data))
colnames(data) <- s@meta.data$merged_cluster
aggregated <- t(rowsum(t(data), group = colnames(data)))
for (col in colnames(aggregated)) {
  aggregated[,col] <- aggregated[,col] / sum(colnames(data) == col)
}

# Subset for marker genes and scale data
plt_data <- aggregated[row.names(aggregated) %in% marker_genes_to_plot,]
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
plt <- plt + geom_dendro(row_clust, pointing = "side", xlim=c(4.5, 6))
plt <- plt + geom_dendro(col_clust, pointing = "updown", ylim=c(19.5, 20))

plt <- plt + labs(x = "Microglia Cluster", y = "Marker Gene") +
  theme(panel.background = element_rect(fill = 'white')) + 
  ggtitle("Scaled Average SCT Expression of FindMarker Genes") + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 2), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_equal() + labs(fill = "Scaled Expression")

# Fig 3M
ggsave(filename = "dendro_heatmap_2_3_4_6.pdf", plot = plt, device = "pdf", path = marker_dir, scale = 0.5, height = 20, width = 15)


