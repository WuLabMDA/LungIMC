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

# celltype_expansion_dir <- file.path(phenotype_dir, "ExpansionInteraction")
# interaction_graph_name <- "expansion_interaction_graph"
# threshold_val <- 50
# cell_type_interaction_name <- paste0("ExpansionInteractionThreshold", threshold_val, ".RData")
# cell_type_interaction_path <- file.path(celltype_expansion_dir, cell_type_interaction_name)
# cell_spatial_path <- file.path(cell_type_interaction_path)
# load(cell_spatial_path)

celltype_delaunay_dir <- file.path(phenotype_dir, "DelaunayInteraction")
interaction_graph_name <- "delaunay_interaction_graph"
threshold_val <- 30
cell_type_interaction_name <- paste0("DelaunayInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_delaunay_dir, cell_type_interaction_name)
load(cell_type_interaction_path)


set.seed(1234)
spe <- aggregateNeighbors(spe, colPairName = interaction_graph_name, 
                          aggregate_by = "metadata", count_by = "celltype")
cn_kmeans <- kmeans(spe$aggregatedNeighbors, centers = 10)
spe$cn_celltypes <- as.factor(cn_kmeans$cluster)

# # save the spe with cellular neighborhood information
# spe_cn_path <- file.path(data_root_dir, "SpatialAnalysis", "spe_cn10.rds")
# saveRDS(spe, spe_cn_path)

# Cellular neighborhood analysis heatmap
for_plot <- prop.table(table(spe$cn_celltypes, spe$celltype), margin = 1)
color_pal <- colorRampPalette(c("dark blue", "white", "dark red"))(10)
pheatmap(for_plot, color = color_pal, scale = "column")