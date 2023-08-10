library(imcRtools)
library(readxl)
library(stringr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(gridExtra)
library(corrplot)
library(RColorBrewer)
library(NCmisc)
library(dittoSeq)
library(viridis)
library(hrbrthemes)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

spe_celltype_name <-"lung_spe_32_cell_subtypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

metadata_dir <- file.path(data_root_dir, "Metadata")
# load ROI information
roi_info_name <- "ROI_Info_Aggregation.csv"
roi_info_path <- file.path(metadata_dir, roi_info_name)
roi_df <- read_csv(roi_info_path)
# load lesion information
slide_info_name <- "Lesion_Info_Aggregation.csv"
slide_info_path <- file.path(metadata_dir, slide_info_name)
slide_info_df <- read_csv(slide_info_path)
# load patient information
patient_info_name <- "Patient_Info.xlsx"
patient_info_path <- file.path(metadata_dir, patient_info_name)
patient_df <- read_excel(patient_info_path)

## filtering ROIs
roi_lst <- roi_df$ROI_ID
interested_roi_vec <- c()
for (ind in 1:nrow(roi_df)) {
    roi_location <- roi_df$ROI_Location[ind]
    roi_name <- roi_df$ROI_ID[ind]
    if (grepl("2571-1D", roi_name))
        next    
    if (roi_location %in% c("Normal", "Tumor")) 
        interested_roi_vec <- append(interested_roi_vec, roi_df$ROI_ID[ind])
}

## obtain cells inside tumor lesions
lesion_indices = c()
cell_id_lst <- rownames(colData(spe))
for (cur_roi in interested_roi_vec) 
    lesion_indices <- append(lesion_indices, which(startsWith(cell_id_lst, cur_roi)))
celltype_lst <- spe$cellsubtype
celltype_lst <- celltype_lst[lesion_indices]
cell_slide_lst <- spe$slide_id[lesion_indices]
cell_stage_lst <- vector("character", length(cell_slide_lst))
for (ind in 1:length(cell_slide_lst)) {
    cell_slide_name <- cell_slide_lst[ind]
    slide_index <- which(slide_info_df$Slide_ID == cell_slide_name)
    cell_stage_lst[ind] <- slide_info_df$Slide_Diag[slide_index]
}

# collect stage cell ratio list
stage_name_lst <- c()
cell_name_lst <- c()
cell_ratio_lst <- c()

# list all cell subtypes
interested_cell_lst <- c("CD163- Macrophages", "Ki67+ Macrophages", "CD163+ Macrophages")
all_stage_lst <- c("Normal", "AAH", "AIS", "MIA", "ADC")


for (cur_stage in all_stage_lst) {
    stage_cell_indices <- which(cell_stage_lst == cur_stage)
    stage_celltype_lst <- celltype_lst[stage_cell_indices]
    stage_celltype_lst <- stage_celltype_lst[stage_celltype_lst %in% interested_cell_lst]
    for (cell_type in interested_cell_lst) {
        cell_ratio <- sum(stage_celltype_lst == cell_type) / length(stage_celltype_lst)
        stage_name_lst <- append(stage_name_lst, cur_stage)
        cell_name_lst <- append(cell_name_lst, cell_type)
        cell_ratio_lst <- append(cell_ratio_lst, cell_ratio)
    }
}

subtype_proportion_dir <- file.path(data_root_dir, "NatureFigures", "Fig02", "ProportionDensity", "ProportionEvolution")
if (!file.exists(subtype_proportion_dir))
    dir.create(subtype_proportion_dir, recursive = TRUE)
stage_subcell_ratio_df <- data.frame(Stage=stage_name_lst, CellType=cell_name_lst, CellRatio=cell_ratio_lst)
stage_subcell_ratio_path <- file.path(subtype_proportion_dir, "stage_macrophage_ratios.csv")
write.csv(stage_subcell_ratio_df, stage_subcell_ratio_path, row.names=FALSE)