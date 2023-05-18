library(SpatialExperiment)
library(imcRtools)
library(cytomapper)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(viridis)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(readxl)
library(knitr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

celltype_delaunay_dir <- file.path(phenotype_dir, "DelaunayInteraction")
interaction_graph_name <- "delaunay_interaction_graph"
threshold_val <- 50
cell_type_interaction_name <- paste0("DelaunayInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_delaunay_dir, cell_type_interaction_name)
load(cell_type_interaction_path)


set.seed(1234)
num_cn <- 8
spe <- aggregateNeighbors(spe, colPairName = interaction_graph_name, 
                          aggregate_by = "metadata", count_by = "celltype")
cn_kmeans <- kmeans(spe$aggregatedNeighbors, centers = num_cn)
spe$cn_celltypes <- as.factor(cn_kmeans$cluster)

# save the spe with cellular neighborhood information
spe_cn_name <- paste0("Delaunay", threshold_val, "-CN", num_cn, ".rds")
spe_cn_path <- file.path(celltype_delaunay_dir, spe_cn_name)
saveRDS(spe, spe_cn_path)

# Cellular neighborhood analysis heatmap
for_plot <- prop.table(table(spe$cn_celltypes, spe$celltype), margin = 1)
color_pal <- colorRampPalette(c("dark blue", "white", "dark red"))(num_cn)
pheatmap(for_plot, color = color_pal, scale = "column")