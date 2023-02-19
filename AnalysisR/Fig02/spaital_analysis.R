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

spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)


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