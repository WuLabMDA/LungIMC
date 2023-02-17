library(imcRtools)
library(SpatialExperiment)
library(stringr)
library(RColorBrewer)
library(dittoSeq)
library(ggplot2)
library(viridis)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

## load SpatialExperiment data
steinbock_spe_name <- "steinbock_spe_2023_02_10"
spe_file_path <- file.path(phenotype_dir, paste0(steinbock_spe_name, ".rds"))
load(spe_file_path)
data$celltype <- data@metadata$celltypes # reorganize cell type information
spe <- data

# steinbock interaction graph
test_roi_name <- "H17-0458-5-ROI016"
plotSpatial(spe[, spe$sample_id == test_roi_name],
            node_color_by = "celltype",
            img_id = "sample_id",
            draw_edges = TRUE,
            colPairName = "neighborhood",
            nodes_first = FALSE,
            edge_color_fix = "grey") +
    ggtitle("steinbock interaction graph")


interaction_out <- testInteractions(spe, group_by = "sample_id",
                                    label = "celltype", 
                                    colPairName = "neighborhood", 
                                    method = "patch",
                                    patch_size = 3,
                                    iter = 200)
cell_type_spatial_path <- file.path(phenotype_dir, paste0("InteractionAnalysisIter200", ".RData"))
save(spe, interaction_out, file = cell_type_spatial_path)