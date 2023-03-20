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

celltype_expansion_dir <- file.path(phenotype_dir, "ExpansionInteraction")
threshold_val <- 24
cell_type_interaction_name <- paste0("ExpansionInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_expansion_dir, cell_type_interaction_name)
cell_spatial_path <- file.path(cell_type_interaction_path)
load(cell_spatial_path)


set.seed(1234)
spe <- aggregateNeighbors(spe, colPairName = "expansion_interaction_graph", 
                          aggregate_by = "metadata", count_by = "celltype")
cn_kmeans <- kmeans(spe$aggregatedNeighbors, centers = 10)
spe$cn_celltypes <- as.factor(cn_kmeans$cluster)

# Cellular neighborhood analysis heatmap
for_plot <- prop.table(table(spe$cn_celltypes, spe$celltype), margin = 1)
color_pal <- colorRampPalette(c("dark blue", "white", "dark red"))(10)
pheatmap(for_plot, color = color_pal, scale = "column")