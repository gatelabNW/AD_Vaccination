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
# Summary: Wilcoxon test to compare nFeatures and MT percent between LCMB and CAA (Ext Fig 3C-D)
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
suppressMessages ({
  library('Seurat')
  library('glue')
  library('dplyr')
  library("tidyverse")
  library('stringr')
  library("Matrix")
  library("data.table")
  library('sctransform')
  library('ggplot2')
  library('glmGamPoi')
})

# Define output folder
output_dir <- "/path/to/general/output/folder/"

# Load integrated Seurat object
load("/path/to/s_integrated_v5")

# Calculate average percent.mt and average nFeatures
meta <- s@meta.data
meta$sample_short <- str_split_fixed(meta$sample_id, "\\.", 3)[,3]
data <- meta %>% dplyr::group_by(sample_short) %>% dplyr::summarize(avg_percent_mt = mean(percent.mt), avg_nFeatures = mean(nFeature_RNA))
data <- dplyr::rename(data, sample = sample_short)
write.csv(data, paste0(output_dir, "full_sample_summary_data.csv", row.names = FALSE))

#-------------------------------------------------------------------------------
# Run Wilcoxon test for MT % and nFeatures

# Load summary data and define broad group variable 
data <- read.csv(paste0(output_dir, "full_sample_summary_data.csv"))
data <- data %>% mutate(group = case_when(str_detect(sample, "A") ~ "A",
                                          str_detect(sample, "B") ~ "B"))

# Calculate Wilcoxon exact p-values (Ext Fig 3C-D)
mt_exact <- stats::wilcox.test(x = data$avg_percent_mt[data$group == "B"], y = data$gw_avg_percent_mt[data$group == "A"],
                               alternative = "two.sided", paired = FALSE, exact = TRUE, correct = FALSE)
nfeat_exact <- stats::wilcox.test(x = data$avg_nFeatures[data$group == "B"], y = data$avg_nFeatures[data$group == "A"],
                                  alternative = "two.sided", paired = FALSE, exact = TRUE, correct = FALSE)

# Calculate Wilcoxon continuity-corrected p-values (not used)
mt_adj <- stats::wilcox.test(x = data$avg_percent_mt[data$group == "B"], y = data$avg_percent_mt[data$group == "A"],
                             alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE)
nfeat_adj <- stats::wilcox.test(x = data$avg_nFeatures[data$group == "B"], y = data$avg_nFeatures[data$group == "A"],
                                alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE)

# Export summary table
summary <- data.frame(value = c("mt_percent", "nfeatures"),
                      w_stat = c(mt_exact$statistic, nfeat_exact$statistic),
                      p_exact = c(mt_exact$p.value, nfeat_exact$p.value),
                      p_corrected = c(mt_adj$p.value, nfeat_adj$p.value))
write.csv(summary, paste0(output_dir, "wilcox_B_vs_A_summary.csv"), row.names = FALSE)

