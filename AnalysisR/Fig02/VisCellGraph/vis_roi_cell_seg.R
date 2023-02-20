library(SpatialExperiment)
library(imcRtools)
library(cytomapper)

set.seed(1234)
# set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
## create spatial analysis directory
spatial_dir <- file.path(data_root_dir, "SpatialAnalysis")

## load spe
spe_name <- "lung_imc_spe"
spe_path <- file.path(spatial_dir, paste0(spe_name, ".rds"))
spe <- readRDS(spe_path)
# load masks
lung_imc_mask_path <- file.path(spatial_dir, "lung_imc_masks.rds")
masks <- readRDS(lung_imc_mask_path)

# Sample images
sample_rois <- sample(unique(spe$sample_id), 2)
sample_masks <- masks[names(masks) %in% sample_rois]

plotCells(sample_masks, object = spe, 
          cell_id = "ObjectNumber", img_id = "sample_id",
          colour_by = "celltype")