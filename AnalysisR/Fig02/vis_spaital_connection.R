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

## load spe object
spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

# test roi
test_roi_name <- "2166-1B-ROI009"

# # cell visualization
# plotSpatial(spe[, spe$sample_id == test_roi_name], 
#             node_color_by = "celltype", 
#             img_id = "sample_id", 
#             node_size_fix = 0.5)

# cell connection visualization
plotSpatial(spe[, spe$sample_id == test_roi_name],
            node_color_by = "celltype",
            img_id = "sample_id",
            node_size_fix = 0.5,
            draw_edges = TRUE,
            colPairName = "neighborhood",
            nodes_first = FALSE,
            directed = FALSE,
            edge_color_fix = "grey")