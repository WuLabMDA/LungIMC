library(imcRtools)
library(ggplot2)
library(ggcorrplot)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
cell_spatial_path <- file.path(phenotype_dir, paste0("InteractionAnalysisIter200", ".RData"))
load(cell_spatial_path)

## load ROI diagnosis information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_path <- file.path(metadata_dir, "ROI_Info.xlsx")
roi_meta_info <- read.xlsx(roi_info_path)


# Normal/AAH/AIS/MIA/ADC
path_stage <- "ADC"
if (path_stage == "Normal") {
    subset_roi_info <- subset(roi_meta_info, ROI_Diag==path_stage)
} else {
    subset_roi_info <- subset(roi_meta_info, ROI_Diag==path_stage & ROI_Location=="Tumor")
}

celltype_order <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell", "CD8-T-Cell", "CD4-T-Cell", 
                    "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte", "MDSC", "Fibroblast", "Undefined")

subset_roi_lst <- subset_roi_info$ROI_ID
subset_roi_num <- length(subset_roi_lst)
cell_id_lst <- rownames(colData(spe))
cell_type_lst <- spe$celltype
roi_cell_ratios <- matrix(nrow=subset_roi_num, ncol=length(celltype_order))


for (ir in 1:subset_roi_num) {
    cur_roi = subset_roi_lst[ir]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- cell_type_lst[cell_indices]
    roi_cell_num <- length(roi_celltypes)
    for (ic in 1:length(celltype_order)) {
        roi_cell_ratios[ir, ic] <- sum(roi_celltypes == celltype_order[ic]) * 1.0 / roi_cell_num
    }
}

cell_ratio_df <- as.data.frame(roi_cell_ratios, row.names = subset_roi_lst)
colnames(cell_ratio_df) <- celltype_order
cell_type_corr <- round(cor(cell_ratio_df), 2)


to_order <- c("Undefined", "Fibroblast", "MDSC", "Monocyte", "Macrophage", "Proliferating-Cell", "T-Reg-Cell", "CD4-T-Cell", "CD8-T-Cell", 
              "Endothelial-Cell", "Dendritic-Cell", "NK-Cell", "Neutrophil", "B-Cell", "Epithelial-Cell")

cell_type_corr <- cell_type_corr[, to_order]
ggcorrplot(cell_type_corr)


