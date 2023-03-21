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
threshold_val <- 30
cell_type_interaction_name <- paste0("DelaunayInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_delaunay_dir, cell_type_interaction_name)
load(cell_type_interaction_path)

# clusting to obtain cn types
set.seed(1234)
spe <- aggregateNeighbors(spe, colPairName = interaction_graph_name, 
                          aggregate_by = "metadata", count_by = "celltype")
cn_kmeans <- kmeans(spe$aggregatedNeighbors, centers = 10)
spe$cn_celltypes <- as.factor(cn_kmeans$cluster)

# save cn
cell_df <- data.frame(cell_id = rownames(spe@colData),
                      cell_type = spe$celltype,
                      cell_cn = spe$cn_celltypes,
                      cell_area = spe$area,
                      cell_maj_ax_len = spe$major_axis_length,
                      cell_min_ax_len = spe$minor_axis_length,
                      cell_eccentricity = spe$eccentricity)
cell_type_path <- file.path(phenotype_dir, "cell_type_cn_morphs.csv")
write.csv(cell_df, cell_type_path, row.names=FALSE)






