library(imcRtools)
library(SpatialExperiment)
library(stringr)
library(RColorBrewer)
library(dittoSeq)
library(ggplot2)
library(viridis)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

celltype_delaunay_dir <- file.path(phenotype_dir, "DelaunayInteraction")
threshold_val <- 50
cell_type_interaction_name <- paste0("DelaunayInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_delaunay_dir, cell_type_interaction_name)
cell_spatial_path <- file.path(cell_type_interaction_path)
load(cell_spatial_path)

# test roi
test_roi_name <- "2166-1B-ROI009"

# cell connection visualization
plotSpatial(spe[, spe$sample_id == test_roi_name],
            node_color_by = "celltype",
            img_id = "sample_id",
            node_size_fix = 0.5,
            draw_edges = TRUE,
            colPairName = "delaunay_interaction_graph",
            nodes_first = FALSE,
            directed = FALSE,
            edge_color_fix = "grey")