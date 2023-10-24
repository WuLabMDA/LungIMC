library(imcRtools)
library(readxl)
library(stringr)
library(tidyverse)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

# cell id
cell_ids <- rownames(colData(spe))
cell_types <- spe$celltype
cell_rois <- colData(spe)$sample_id

# load ROI information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_name <- "ROI_Info_Aggregation.csv"
roi_info_path <- file.path(metadata_dir, roi_info_name)
roi_df <- read_csv(roi_info_path)
cell_stages <- vector("character", length(cell_ids))
for (ind in 1:length(cell_rois)) {
    cell_roi <- cell_rois[ind]
    roi_index <- which(roi_df$ROI_ID == cell_roi)
    cell_stages[ind] <- roi_df$ROI_Diag[roi_index]
}

TIM3_vals <- as.numeric(counts(spe)["TIM3",])
ICOS_vals <- as.numeric(counts(spe)["ICOS",])
TIGIT_vals <- as.numeric(counts(spe)["TIGIT",])
CD73_vals <- as.numeric(counts(spe)["CD73",])
PDL1_vals <- as.numeric(counts(spe)["PDL1",])
LAG3_vals <- as.numeric(counts(spe)["LAG3",])
VISTA_vals <- as.numeric(counts(spe)["VISTA",])
PD1_vals <- as.numeric(counts(spe)["PD1",])
B7H3_vals <- as.numeric(counts(spe)["B7H3",])
CTLA4_vals <- as.numeric(counts(spe)["CTLA4",])

# contruct dataset
cell_exp_df <- data.frame(cell_id = cell_ids, cell_roi = cell_rois, 
                          cell_type = cell_types, cell_stage = cell_stages, 
                          TIM3 = TIM3_vals, ICOS = ICOS_vals, TIGIT = TIGIT_vals, 
                          CD73 = CD73_vals, PDL1 = PDL1_vals, LAG3 = LAG3_vals, 
                          VISTA = VISTA_vals, PD1 = PD1_vals,
                          B7H3 = B7H3_vals, CTLA4 = CTLA4_vals)

# save data
tim3_dir <- file.path(phenotype_dir, "TIM3")
if (!file.exists(tim3_dir))
    dir.create(tim3_dir, recursive = TRUE)
tim3_path <- file.path(tim3_dir, "cell_tim3.rds")
saveRDS(cell_exp_df, file = tim3_path)




