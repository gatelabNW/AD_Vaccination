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
# Summary: Cohort 1 gene correlations with distance to nearest amyloid plaque
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
  library(ggpubr)
})

options(future.globals.maxSize=1048576000000)
load("helperfunctions.RData")
setwd("/path/to/data")

# Load seurat object
s <- readRDS("/path/to/seurat.rds")

# Subset to gray matter
s <- subset(s, manual_annotation %!in% c("meninges", "white"))

# Subset for amyloid filter
s <- subset(s, amyloid_filter == "include")

# Flag NA amy spots
s@meta.data$amyinclude <- ifelse(is.na(s@meta.data$amyloid_fluo), "Exclude", "Include")

# Subset to non-NA spots
s <- subset(s, amyinclude == "Include")

# Lists to keep track of total amyloid per sample and samples to remove
countlist <- samplestoremove <- list()

# Loop over sample IDs
for(sample in unique(s@meta.data$sample_id)){
  
  # Sum all amyloid for this sample
  countlist[[sample]] <- sum(s@meta.data$amyloid_fluo[s@meta.data$sample_id == sample])
  
  # If amy sum = 0, store sample ID in samplestoremove
  samplestoremove[[sample]] <- ifelse(countlist[[sample]] == 0, sample, 0)
  
}

# rbind samplestoremove
srm <- list.rbind(samplestoremove)

# filter to IDs to be removed
srm <- srm[srm != 0]

# subset to samples with amyloid > 0
s <- subset(s, sample_id %!in% srm)

# Initialize sample meta list
smeta <- list()

# loop over sample IDs
for(name in unique(s@meta.data$sample_id)){
  
  # Take sample-specific meta data
  smeta[[name]] <- s@meta.data[s@meta.data$sample_id == name,]
  
  # Add row
  smeta[[name]]$row <- s@images[[name]]@coordinates$row
  
  # Add col
  smeta[[name]]$col <- s@images[[name]]@coordinates$col
  
  # Add binary var for amyloid
  smeta[[name]]$amyspot <- ifelse(smeta[[name]]$amyloid_fluo > 125, "Amyloid", "Not Amyloid")
  
  # Add region var
  smeta[[name]]$region <- ifelse(smeta[[name]]$manual_annotation == "meninges", "Meninges", "Gray/White")
  
}

# Loop over sample IDs
for(sample in names(smeta)){
  
  # store smeta[[sample]] in dat
  dat <- smeta[[sample]]
  
  # initialize distance var
  dat$distance <- smeta[[sample]]$distance <- 0
  
  # loop over all spots that are not amyloid
  for(spot in unique(dat$barcode[dat$amyspot != "Amyloid"])){
    
    # current row
    crow <- dat$row[dat$barcode == spot]
    
    # current column
    ccol <- dat$col[dat$barcode == spot]
    
    # current region
    clayer <- dat$region[dat$barcode == spot]
    
    # modified distance formula to accomodate Seurat row/col convention
    dat$distance <- sqrt((dat$row - crow)^2 + ((dat$col - ccol)/2)^2)
    
    if(sample != "A18.148"){
      # distance = minimum distance to amyloid plaque in the same region
      smeta[[sample]]$distance[smeta[[sample]]$barcode == spot] <- min(dat$distance[dat$amyspot == "Amyloid" & dat$region == clayer])
    }else{
      # A18-148 has no amyspots in meninges and returns infinity - relax layer restriction for this sample only
      smeta[[sample]]$distance[smeta[[sample]]$barcode == spot] <- min(dat$distance[dat$amyspot == "Amyloid"])
    }
    
  }
  
  # store max distance
  mdist <- max(smeta[[sample]]$distance)
  
  # invert distances per sample - turn max into zero and min into new max
  smeta[[sample]]$Distance_Inverted <- mdist - smeta[[sample]]$distance
  
  # calculate scaled distance
  smeta[[sample]]$distance_scaled <- (smeta[[sample]]$distance - min(smeta[[sample]]$distance)) / (max(smeta[[sample]]$distance) - min(smeta[[sample]]$distance))*100
  
  # calculate scaled inverted distance
  smeta[[sample]]$distance_inv_scaled <- (smeta[[sample]]$Distance_Inverted - min(smeta[[sample]]$Distance_Inverted)) / (max(smeta[[sample]]$Distance_Inverted) - min(smeta[[sample]]$Distance_Inverted))*100
  
  # calculate scaled amyloid
  smeta[[sample]]$amyloid_fluo_scaled <- (smeta[[sample]]$amyloid_fluo - min(smeta[[sample]]$amyloid_fluo)) / (max(smeta[[sample]]$amyloid_fluo) - min(smeta[[sample]]$amyloid_fluo))*100
  
  
}

# Add var for C2L microglia data
s@meta.data$microglia <- 0

# Initialize list for microglia info
gliadat <- list()

# copy smeta
smeta2 <- smeta

# loop over samples to attach microglia data
for(sample in names(smeta2)){
  
  # read csv
  gliadat[[sample]] <- read.csv(paste0("/output_dir/" , sample, "_combined_meta.csv"))
  
  # remove NAs
  gliadat[[sample]] <- gliadat[[sample]][!is.na(gliadat[[sample]]$Microglia),]
  
  # ensure same cells
  gliadat[[sample]] <- gliadat[[sample]][gliadat[[sample]]$barcode %in% smeta2[[sample]]$barcode,]
  
  # order by barcode
  gliadat[[sample]] <- gliadat[[sample]][order(gliadat[[sample]]$barcode),]
  
  # ensure same cells
  smeta2[[sample]] <- smeta2[[sample]][smeta2[[sample]]$barcode %in% gliadat[[sample]]$barcode,]
  
  # order by barcode
  smeta2[[sample]] <- smeta2[[sample]][order(smeta2[[sample]]$barcode),]
  
  # ensure all barcodes match
  print(all.equal(gliadat[[sample]]$barcode, smeta2[[sample]]$barcode))
  
  # attach microglia data
  smeta2[[sample]]$microglia <- gliadat[[sample]]$Microglia
  
  # add scaled microglia data
  smeta2[[sample]]$microglia_scaled <- (smeta2[[sample]]$microglia - min(smeta2[[sample]]$microglia)) / (max(smeta2[[sample]]$microglia) - min(smeta2[[sample]]$microglia))*100
  
}

# rbind smeta
newmeta <- list.rbind(smeta2)

# store unchanged meta data
oldmeta <- s@meta.data

# subset to spots with microglia data
s <- subset(s, sample_barcode %in% newmeta$sample_barcode)

# meta = newmeta
s@meta.data <- newmeta

# Set rownames
rownames(s@meta.data) <- s@meta.data$sample_barcode

# Set default assay to SCT
DefaultAssay(s) <- "SCT"

# make sure only nAD, lim and ext are included
s <- subset(s, sample_id %in% names(smeta))

# Make SCT model assays all the same
for (name in names(s@assays$SCT@SCTModel.list)){
  
  s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
  
}

# Prep for FindMarkers
s <- PrepSCTFindMarkers(object = s)

# Set condition_clearance as Idents for upcoming FoldChange()s
Idents(s) <- s@meta.data$condition_clearance

# Set ident 1 and 2 to comparison groups
LFC1 <- FoldChange(object = s, ident.1 = "nAD", ident.2 = "lim", fc.name = "avg_log2FC", base = 2) %>% data.frame()

# Keep genes expressed in at least 1% of either group
genes_keep1 <- row.names(LFC1)[LFC1$pct.1 >= 0.01 | LFC1$pct.2 >= 0.01]

# Set ident 1 and 2 to comparison groups
LFC2 <- FoldChange(object = s, ident.1 = "nAD", ident.2 = "ext", fc.name = "avg_log2FC", base = 2) %>% data.frame()

# Keep genes expressed in at least 1% of either group
genes_keep2 <- row.names(LFC2)[LFC2$pct.1 >= 0.01 | LFC2$pct.2 >= 0.01]

# Combine lists without duplicates
genes_keep <- unique(c(genes_keep1, genes_keep2))

# Calculate CDR
cdr <- colMeans(GetAssayData(s, assay = "SCT", layer = "data") > 0)

# Add CDR to meta
s@meta.data$cdr <- cdr

# Add centered cdr to meta
s@meta.data$cdr_centered <- (s@meta.data$cdr - mean(s@meta.data$cdr)) / sd(s@meta.data$cdr)

# Add centered age to meta
s@meta.data$age_centered <- (s@meta.data$age - mean(s@meta.data$age)) / sd(s@meta.data$age)

# Add centered gDNA to meta
s@meta.data$gDNA_percent_centered <- (s@meta.data$gDNA_percent - mean(s@meta.data$gDNA_percent)) / sd(s@meta.data$gDNA_percent)

# Initalize lists for metadata and assay data
smeta <- dat <- list()

# Get whole assay data
alldat <- as.data.frame(t(GetAssayData(s, assay = "SCT", layer = "data")))

# These genes are ribosomal, mitochondrial, or hemoglobin
rcols <- grep(pattern = "^RPS|^RPL|^MT-|^HB", x = colnames(alldat))

# Remove RP, MT and HB genes
alldat <- alldat[,-c(rcols)]

# Remove and ribosomal, mitochondrial, or hemoglobin genes from genes_keep
genes_keep <- genes_keep[genes_keep %in% colnames(alldat)]

# Filter to genes_keep
alldat <- alldat[, genes_keep]

# Make sure all spots are accounted for
sum(rownames(alldat) %in% s@meta.data$sample_barcode) == nrow(s@meta.data)

# Populate meta and assay lists
for(name in unique(s@meta.data$sample_id)){
  
  smeta[[name]] <- s@meta.data[s@meta.data$sample_id == name,]
  
  smeta[[name]] <- smeta[[name]][order(smeta[[name]]$sample_barcode),]
  
  dat[[name]] <- alldat[rownames(alldat) %in% smeta[[name]]$sample_barcode,]
  
}

# free unused memory
gc()


# Use all genes
genes <- colnames(alldat)

# Populate mats list by group
mats <- list(
  
  "nAD" = rbind(dat[["A34583.2"]][,genes], dat[["N35127N"]][,genes], dat[["A34995"]][,genes],
                dat[["A34993.2"]][,genes]),
  
  "Limited" = rbind(dat[["AN1792.102.10"]][,genes], dat[["AN1792.102.11"]][,genes], dat[["AN1792.102.17"]][,genes],
                    dat[["AN1792.102.19"]][,genes], dat[["AN1792.102.22"]][,genes], dat[["AN1792.102.4"]][,genes]),
  
  "Extensive" = rbind(dat[["AN1792.102.1"]][,genes], dat[["AN1792.102.16"]][,genes], dat[["AN1792.102.20"]][,genes], dat[["AN1792.102.21"]][,genes],
                      dat[["AN1792.102.6"]][,genes], dat[["AN1792.102.7"]][,genes], dat[["AN1792.102.8"]][,genes])
  
)

# Populate metas list by group
metas <- list(
  
  "nAD" = rbind(smeta[["A34583.2"]], smeta[["N35127N"]], smeta[["A34995"]],
                smeta[["A34993.2"]]), 
  
  "Limited" = rbind(smeta[["AN1792.102.10"]], smeta[["AN1792.102.11"]], smeta[["AN1792.102.17"]],
                    smeta[["AN1792.102.19"]], smeta[["AN1792.102.22"]], smeta[["AN1792.102.4"]]),
  
  "Extensive" = rbind(smeta[["AN1792.102.1"]], smeta[["AN1792.102.16"]], smeta[["AN1792.102.20"]], smeta[["AN1792.102.21"]],
                      smeta[["AN1792.102.6"]], smeta[["AN1792.102.7"]], smeta[["AN1792.102.8"]])
  
)

# Initialize lists to be used in correlations
resmat <- geneset <- covs <- list()

# Free unused memory
gc()

# Loop over conditions to use the lists we just made
for(condition_group in names(mats)){
  
  # make a results DF that will adapt its size to fit the data we're analyzing
  mat <- data.frame(t(mats[[condition_group]][1:2,]))
  
  # We only want the shape, not the data
  mat[,c(1:2)] <- 0
  
  # set colnames
  colnames(mat) <- c("Spearmans_Rho", "p")
  
  # Covariate matrix for partial Spearman cor tests
  covs[[condition_group]] <- metas[[condition_group]][,"cdr_centered"]
  
  # print group start time
  print(paste0(condition_group, " start: ", Sys.time()))
  
  # Now loop over each gene and calculate covariate-adjusted Spearman cor
  for(gene in rownames(mat)){
    
    # Partial Spearman cors
    ctest <- pcor.test(x = metas[[condition_group]]$distance, 
                       y = mats[[paste0(condition_group)]][,paste0(gene)], 
                       z = covs[[condition_group]],
                       method = "spearman")
    
    # Put rho in col 1
    mat[paste0(gene),"Spearmans_Rho"] <- ctest$estimate
    
    # Put pval in col 2
    mat[paste0(gene),"p"] <- ctest$p.value
    
  }
  
  # Print group end time
  print(paste0(condition_group, " end: ", Sys.time()))
  
  # make sure genes are rownames
  mat$gene <- rownames(mat)
  
  # should be zero
  print(paste0("NA sum for group ", condition_group, ": ", sum(is.na(mat))))
  
  # calculate BH adjusted pvals
  mat$padj <- p.adjust(mat$p, method="BH")
  
  # order mat by smallest to largest padj
  mat <- data.frame(mat[order(mat$padj),])
  
  # remove any NAs
  mat <- na.omit(mat)
  
  # binary var for padj sig
  mat$padjsig <- ifelse(mat$padj < 0.05, "Sig", "Not Sig")
  
  # binary var for rho sig
  mat$rhosig <- ifelse(abs(mat$Spearmans_Rho) > 0.1, "Sig", "Not Sig")
  
  # mat is now resmat
  resmat[[paste0(condition_group)]] <- mat
  
  # Add a group var for upset plots
  resmat[[paste0(condition_group)]]$Group <- condition_group
  
  # This may be redundant
  resmat[[paste0(condition_group)]] <- resmat[[paste0(condition_group)]][order(resmat[[paste0(condition_group)]]$padj),]
  
  # Save CSV of cors
  write.csv(resmat[[paste0(condition_group)]], file = paste0(condition_group, "_Amyloid_Corr.csv"))
  
  # Store geneset for each group in a different list
  geneset[[paste0(condition_group)]] <- resmat[[paste0(condition_group)]]
  
}

# Separate lists for labeling
geneset2 <- geneset3 <- geneset4 <- posgenes <- neggenes <- list()

# loop over list names
for(name in names(geneset)){
  
  # geneset 2 is all sig genes
  geneset2[[name]] <- geneset[[name]][geneset[[name]]$padjsig == "Sig" & geneset[[name]]$rhosig == "Sig",]
  
  # geneset 3 is all sig genes, only names - for upset plot
  geneset3[[name]] <- geneset2[[name]]$gene
  
  # geneset 4 is all non-sig genes
  geneset4[[name]] <- geneset[[name]][geneset[[name]]$padjsig != "Sig" | geneset[[name]]$rhosig != "Sig",]
  
  # posgenes is all positively correlated sig genes
  posgenes[[name]] <- geneset2[[name]][geneset2[[name]]$Spearmans_Rho > 0,]
  
  # neggenes is all negatively correlated sig genes
  neggenes[[name]] <- geneset2[[name]][geneset2[[name]]$Spearmans_Rho < 0,]
  
}

# sig genes ls for upset
sig_genes_ls <- geneset3

# only use groups with more than 0 sig genes
sig_genes_ls <- sig_genes_ls[lapply(sig_genes_ls, length) > 0]

# summary DF
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Set color palette for upset plot
suppressWarnings(
  colrs <- rep(c(RColorBrewer::brewer.pal(n = length(sig_genes_ls), name = "Dark2")), 3))

# set color palette
bar_colors <- distinctColorPalette(7)

# initialize percent_DE
percent_DE <- c()

# clust_results is a list containing DE results for each cluster for a given comparison (ex: iAD vs. nAD)
for (cluster in names(sig_genes_ls)) {
  percent <- length(sig_genes_ls[[cluster]]) / ncol(mats[[cluster]])
  percent <- round(percent*100, 2)
  percent_DE <- append(percent_DE, paste0(cluster, ": ", percent, "%"))
}

# clust_genes is a list containing DEGs for each cluster for a given comparison (ex: iAD vs. nAD)
names(sig_genes_ls) <- percent_DE 

# Make upset plot
p <- upset(fromList(sig_genes_ls),
           nsets = length(sig_genes_ls), set_size.show = FALSE, text.scale = 1.5, 
           order.by = "freq", shade.color = "black", point.size = 3, line.size = 1, nintersects = 20,
           sets.bar.color = "darkblue", empty.intersections = "off", main.bar.color = bar_colors)

# Print upset plot
pdf(file = "AmyloidSigGenesUpset.pdf", width = 12, height = 8)
print(p)
dev.off()

# load gene sets
eggen <- read.csv("eggen-microglia.csv")
eggen_genes <- eggen$Gene
ham <- read.csv("hamgenes.csv")
ham_genes <- ham$gene
additional_genes <- c("APOE", "TREM2", "SPP1", "C1QB", "TYROBP", "CD74")
combined_genes <- unique(c(eggen_genes, ham_genes))

# sig genes
newset <- geneset2

# not sig
nset <- geneset4

for(comp in names(newset)){
  
  # Initialize gene label var
  newset[[comp]]$genelabel <- NA
  
  # Specify correlation direction
  newset[[comp]]$direction <- ifelse(newset[[comp]]$Spearmans_Rho >= 0, "Positive", "Negative")
  
  # Specify correlation direction
  nset[[comp]]$direction <- ifelse(nset[[comp]]$Spearmans_Rho >= 0, "Positive", "Negative")
  
  # Calculate PFC
  newset[[comp]]$PFC <- -log10(newset[[comp]]$padj) * abs(newset[[comp]]$Spearmans_Rho)
  
  # Order by PFC
  newset[[comp]] <- newset[[comp]][order(-newset[[comp]]$PFC),]
  
  # Label top 10 PFC genes
  newset[[comp]]$genelabel[1:10] <- newset[[comp]]$gene[1:10]
  
  # hold them to the side
  sigset <- newset[[comp]][1:10,]
  
  # tempset is for labeling from combined and additional gene lists
  tempset <- newset[[comp]][11:nrow(newset[[comp]]),]
  
  # Set ngenes in case there are less than 20 genes
  ngenes <- ifelse(nrow(tempset[tempset$gene %in% combined_genes,]) < 20, nrow(tempset[tempset$gene %in% combined_genes,]), 20)
  
  # add combined genes labels
  tempset$genelabel[tempset$gene %in% combined_genes][1:ngenes] <- tempset$gene[tempset$gene %in% combined_genes][1:ngenes]
  
  # stitch back together
  newset[[comp]] <- rbind(sigset, tempset)
  
  # Force labeling from additional genes
  newset[[comp]]$genelabel[newset[[comp]]$gene %in% additional_genes] <- newset[[comp]]$gene[newset[[comp]]$gene %in% additional_genes]
  
  # nset isn't sig so no label
  nset[[comp]]$genelabel <- NA
  
  # Calculate PFC
  nset[[comp]]$PFC <- -log10(nset[[comp]]$padj) * abs(nset[[comp]]$Spearmans_Rho)
  
  # Only sig if both padj and rho are also sig
  newset[[comp]]$sig <- ifelse(newset[[comp]]$rhosig == "Sig" & newset[[comp]]$padjsig == "Sig", "Sig", "Not Sig")
  
  # Not sig
  nset[[comp]]$sig <- "Not Sig"
  
  # sharedUnique will designate whether a gene correlation is also sig in the same direction in nAD
  nset[[comp]]$sharedUnique <- NA
  
  # nAD gets NA
  if(comp == "nAD"){
    
    newset[[comp]]$sharedUnique <- NA
    
  }else{
    
    # For positive corrs
    newset[[comp]]$sharedUnique[newset[[comp]]$direction == "Positive"] <- ifelse(newset[[comp]]$gene[newset[[comp]]$direction == "Positive"] %in% posgenes[["nAD"]]$gene, "Sig in nAD", NA)
    
    # For negative corrs
    newset[[comp]]$sharedUnique[newset[[comp]]$direction == "Negative"] <- ifelse(newset[[comp]]$gene[newset[[comp]]$direction == "Negative"] %in% neggenes[["nAD"]]$gene, "Sig in nAD", NA)
    
  }
  
  # rbind the two sets so there is one set for the comp
  newset[[comp]] <- rbind(newset[[comp]], nset[[comp]])
  
}

# rbind all comp sets together to form plotting DF
plotgenes <- list.rbind(newset)

# - log(padj) for plotting
plotgenes$logpadj <- -log10(plotgenes$padj)

# Order descending by PFC
plotgenes <- plotgenes[order(-plotgenes$PFC),]

# If the label isn't "sig in nAD", it needs to be the group. For plot faceting 
plotgenes$sharedUnique[is.na(plotgenes$sharedUnique)] <- plotgenes$Group[is.na(plotgenes$sharedUnique)]

# remove nAD
plotgenes <- plotgenes[plotgenes$Group != "nAD",]

# Create new var for plot coloring
plotgenes$plotcol <- plotgenes$sharedUnique

# We're coloring by significance, not sample
plotgenes$plotcol[plotgenes$sig == "Not Sig"] <- "Not Sig"

# for faceting
plotgenes$sharedUnique[plotgenes$sharedUnique == "Sig in nAD"] <- plotgenes$Group[plotgenes$sharedUnique == "Sig in nAD"]

# for coloring non-sig spots gray
plotgenes$sharedUnique[plotgenes$sig == "Not Sig"] <- "Not Sig"

# Only do this if there are very few sig genes and the labels won't become overwhelming
plotgenes$genelabel[plotgenes$sig == "Sig"] <- plotgenes$gene[plotgenes$sig == "Sig"]

# for labeling positive and negative correlations
plotgenes$sharedUnique2 <- plotgenes$sharedUnique

# Label positive correlations
plotgenes$sharedUnique2[plotgenes$Spearmans_Rho < 0 & plotgenes$sig == "Sig"] <- "Positively Correlated"

# Label negative correlations
plotgenes$sharedUnique2[plotgenes$Spearmans_Rho > 0 & plotgenes$sig == "Sig"] <- "Negatively Correlated"

# Set color vector
colvec2 <- c("Positively Correlated" = "red",
             "Negatively Correlated" = "blue",
             "Sig in nAD" = "black",
             "Not Sig" = "gray")

plotgenes$sharedUnique3 <- plotgenes$sharedUnique2

plotgenes$sharedUnique3[plotgenes$plotcol == "Sig in nAD"] <- "Sig in nAD"

plotgenes$sharedUnique3[is.na(plotgenes$sharedUnique3)] <- plotgenes$sharedUnique2[is.na(plotgenes$sharedUnique3)]

# Make plot
p <- ggplot(data = plotgenes, 
            aes(x = -Spearmans_Rho, 
                y = logpadj, 
                label = genelabel)) + 
  geom_point(data = plotgenes, 
             aes(color = sharedUnique3, 
                 size = PFC), alpha = 0.75) + 
  theme_cowplot() +
  geom_label_repel(data = plotgenes, 
                   fill = alpha(c("white"),0.01),  
                   min.segment.length = 0.1,
                   aes(segment.alpha = 0.6, 
                       color = sharedUnique3),
                   force = 15,
                   force_pull = 1,
                   max.overlaps = Inf,
                   max.iter = 50000,
                   max.time = 5,
                   direction = "both",
                   label.padding = 0.2,
                   size = 8,
                   fontface = "bold",
                   show.legend = FALSE,
                   label.size = NA
  ) +
  theme(
    plot.background = element_rect(fill = "white"),
    axis.line.x.bottom=element_line(color="black"),
    axis.line.y.left=element_line(color="black"),
    axis.title = element_text(color = "black", size = 20),
    legend.text = element_text(color="black", size = 16),
    plot.title = element_text(color="black", size = 16),
    title = element_text(color="black"),
    axis.ticks = element_line(color="black"),
    axis.text = element_text(color="black", size = 16)
  ) +
  facet_wrap(.~factor(Group, levels = c("Limited", "Extensive")), scales = "free_y") +
  ylab("-log10(padj)") + 
  xlab("- Spearman's Rho") + 
  ggtitle("All Aß Genes | SCT Data | Distance from Nearest Plaque | Gray Matter | 125") +
  scale_color_manual(values = colvec2) +
  geom_hline(yintercept = 1.30103, 
             color = "lightgray", 
             linetype = "dashed") +
  geom_vline(xintercept = 0.1, 
             color = "lightgray", 
             linetype = "dashed") +
  geom_vline(xintercept = -0.1, 
             color = "lightgray", 
             linetype = "dashed")

pdf(file = "AmyloidSigGenesCombined.pdf", width = 24, height = 12)
print(p)
dev.off()

write.csv(plotgenes, file = "AmyloidSigGenes.csv")
