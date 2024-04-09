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
# Summary: UpSet plot for DESeq2 and MAST results
#
#-------------------------------------------------------------------------------

# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggrepel")
  library("UpSetR")
  library("randomcoloR")
  library("stringr")
})


# Define path to folder containing csv files of all results to use in upset plot
# Ensure that the only files in this folder are the results, and the names correspond to desired set names in plot
result_dir <- "/path/to/all/results/folder/"

# Set these variables to column names corresponding to LFC, adjusted p-value and gene name (so they can be renamed in the plot function)
lfc <- "<name of LFC column>"
padj <- "<name of padj column>"
gene <- "<name of gene name column>"

# Define DE thresholds
fc_thresh <- log2(1.5) 
p_thresh <- 0.05

#-------------------------------------------------------------------------------
# Function to get list of DEGs to use for UpSet plot
load_degs <- function(result_dir, pct_in_title = FALSE, n_in_title = FALSE) {
  
  # Initialize list for DEG sets 
  deg_list <- list()
  
  # Load all results and add DEGs to list
  for (file in list.files(result_dir)) {
    
    # Load results
    results <- read.csv(paste0(result_dir, file), row.names = 1)
    
    # Rename variables 
    results$lfc <- results[[lfc]]
    results$padj <- results[[padj]]
    results$gene <- results[[gene]]
    
    # Define DEGs
    results$DE <- "Not DE"
    results$DE[results$padj < p_thresh & abs(results$lfc) > fc_thresh] <- "DE"
    
    # Add to upset list
    set_name <- str_replace_all(file, ".csv", "")
    degs <- results$gene[results$DE == "DE"]
    
    # Customize set name
    if (pct_in_title & n_in_title) {
      stop("Specify only one metric to include in title!")
    } else if (pct_in_title) {
      pct <- round((length(degs)/nrow(results))*100, 2)
      pct <- paste0(pct, " %")
      deg_list[[paste0(set_name, ": ", pct, " DE")]] <- degs
    } else if (n_in_title) {
      deg_list[[paste0(set_name, ": ", length(degs), " DEGs")]] <- degs
    } else {
      deg_list[[set_name]] <- degs
    }
  }
  
  return(deg_list)
}

# Generate set list for upset with basic set names (set pct_in_title or n_in_title to TRUE to include percent DE or number of DEGs in set names)
deg_list <- load_degs(result_dir)

# Generate initial plot to extract number of intersection bars (needed to generate colors)
plt <- UpSetR::upset(UpSetR::fromList(deg_list), nsets = length(deg_list))
n_bars <- nrow(unique(plt$New_data))

# Adjust seed until colors look good (or use another method to generate n_bars colors)
set.seed(125)
bar_colors <- randomcoloR::distinctColorPalette(n_bars)

# Generate UpSet plot (can customize sets.bar.color, nintersects, etc)
plt <- UpSetR::upset(fromList(deg_list), order.by = "freq", nsets = length(deg_list), nintersects = 35,
                     sets.bar.color = "darkblue", main.bar.color = bar_colors)

