library(imcRtools)
library(tidyverse)
library(readxl)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

## save updated data
updated_spe_name <- "lung_spe_32_cell_subtypes_final"
updated_spe_path <- file.path(phenotype_dir, paste0(updated_spe_name, ".rds"))

# obtain cell
spe <- readRDS(updated_spe_path)
cell_roi_cell_ids <- spe$ObjectNumber
cell_celltypes <- spe$celltype
cell_cellsubtypes <- spe$cellsubtype

# prepare directory
celltype_dir <- file.path(phenotype_dir, "ROI-CellTypes")
if (!dir.exists(celltype_dir))
    dir.create(celltype_dir)

# traverse each roi
cell_rois <- spe$sample_id
unique_roi_lst <- unique(cell_rois)
for (cur_roi in unique_roi_lst) {
    roi_cell_indices <- which(cell_rois == cur_roi) 
    roi_cell_ids <- cell_roi_cell_ids[roi_cell_indices]
    roi_celltypes <- cell_celltypes[roi_cell_indices]
    roi_cellsubtypes <- cell_cellsubtypes[roi_cell_indices]

    cur_roi_df <- data.frame(cell_id = roi_cell_ids,
                             celltype = roi_celltypes,
                             cellsubtype = roi_cellsubtypes)
    cur_roi_path <- file.path(celltype_dir, paste0(cur_roi, ".csv"))
    write.csv(cur_roi_df, cur_roi_path, row.names=FALSE)
}


