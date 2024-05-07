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
# Written by: Anne Forsyth  & Thomas Watson
# Summary: Dot plot for pathways enriched in microglia-enriched gray matter (Fig 2F, Ext Fig 2A)
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("enrichR")
  library("clusterProfiler")
  library("pathview")
  library("ggplot2")
  library("rlist")
})

# Define output folder
output_folder <- "/path/to/general/gseKEGG/output/folder/"

# Load gseKEGG results
lim_nAD <- readRDS(paste0(output_folder, "gseKEGG_all_genes/microglia_gray/lim_nAD_adj.rds"))
ext_nAD <- readRDS(paste0(output_folder, "gseKEGG_all_genes/microglia_gray/ext_nAD_adj.rds"))

#-------------------------------------------------------------------------------
# Format data and generate dot plot

# Extract results from RDS objects 
lim <- lim_nAD@result
ext <- ext_nAD@result

# Create variable for comparison name
lim$comp <- "iAD-Lim vs. nAD"
ext$comp <- "iAD-Ext vs. nAD"

# Combine comparison data
df <- rbind(lim, ext)

# Define variable for dot color
df$plotcol <- 0
df$plotcol[df$comp == "iAD-Lim vs. nAD"] <- "olivedrab3"
df$plotcol[df$comp == "iAD-Ext vs. nAD"] <- "cadetblue3"

# Order results by enrichment score
df$plotorder <- rank(df$enrichmentScore)
df <- df[order(df$plotorder),]

# Generate dot plot
p1 <- ggplot(df, aes(x = enrichmentScore, y = reorder(Description, plotorder), fill = comp, label = Description)) +
  geom_point(size = 5, stat="identity", stroke = NA, shape = 21) + 
  geom_text(size = 4.25, vjust = -100) +
  labs(fill = "Comparison") + ggtitle("gseKEGG Analysis: Microglia-Gray") +
  theme_bw() + ylab("") + xlab("Enrichment Score") +
  scale_fill_manual(values = c("olivedrab3", "cadetblue3")) +
  theme(axis.text.y = element_text(color = df$plotcol)) + xlim(c(-.9,.9)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted")

# Fig 2F, Ext Fig 2A
pdf(paste0(output_folder, "gseKEGG_all_genes/microglia_gray_final/lim_ext_plt.pdf"), width = 8, height = 10)
print(p1)
dev.off()

