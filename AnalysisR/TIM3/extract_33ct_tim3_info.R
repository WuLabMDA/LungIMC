library(readxl)
library(stringr)
library(tidyverse)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

spe_celltype_name <-"lung_spe_32_cell_subtypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

# cell id
cell_ids <- rownames(colData(spe))
cell_types <- spe$cellsubtype
cell_rois <- colData(spe)$sample_id

# load ROI information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_name <- "ROI_Info_Aggregation.csv"
roi_info_path <- file.path(metadata_dir, roi_info_name)
roi_df <- read_csv(roi_info_path)
cell_stages <- vector("character", length(cell_ids))
cell_locations <- vector("character", length(cell_ids))
for (ind in 1:length(cell_rois)) {
  cell_roi <- cell_rois[ind]
  roi_index <- which(roi_df$ROI_ID == cell_roi)
  cell_stages[ind] <- roi_df$ROI_Diag[roi_index]
  cell_locations[ind] <- roi_df$ROI_Location[roi_index]
}

TIM3_vals <- as.numeric(counts(spe)["TIM3",])

# contruct dataset
cell_exp_df <- data.frame(cell_id = cell_ids, cell_roi = cell_rois, cell_type = cell_types, 
                          cell_location = cell_locations, cell_stage = cell_stages, TIM3 = TIM3_vals)
cell_exp_df <- cell_exp_df[cell_exp_df$cell_location != "AdjacentNormal",]

# save data
tim3_dir <- file.path(phenotype_dir, "TIM3")
if (!file.exists(tim3_dir))
  dir.create(tim3_dir, recursive = TRUE)
tim3_path <- file.path(tim3_dir, "cell_subtype_tim3.rds")
saveRDS(cell_exp_df, file = tim3_path)
