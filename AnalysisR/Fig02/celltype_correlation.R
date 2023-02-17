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


celltype_order <- c("Epithelial", "B cells", "Neutrophils", "NK cell", "Dendritic cell", 
                    "Endothelial cells", "CD8 T cells", "CD4 T cell", "T-reg cells", "Macrophage", 
                    "Monocytes", "MDSC", "Smooth muscle/Stromal", "Unknown")
## find all unique sample ids
unique_roi_lst <-unique(spe$sample_id)
uniuqe_roi_num <- length(unique_roi_lst)
cell_id_lst <- rownames(colData(spe))
cell_type_lst <- spe$celltype
roi_cell_ratios <- matrix(nrow=uniuqe_roi_num, ncol=length(celltype_order))

ttl_cell_num <- 0
for (ir in 1:uniuqe_roi_num) {
    cur_roi = unique_roi_lst[ir]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- cell_type_lst[cell_indices]
    roi_cell_num <- length(roi_celltypes)
    ttl_cell_num <- ttl_cell_num + roi_cell_num
    for (ic in 1:length(celltype_order)) {
        roi_cell_ratios[ir, ic] <- sum(roi_celltypes == celltype_order[ic]) * 1.0 / roi_cell_num
    }
}

cell_ratio_df <- as.data.frame(roi_cell_ratios, row.names = unique_roi_lst)
colnames(cell_ratio_df) <- celltype_order
cell_type_corr <- round(cor(cell_ratio_df), 2)
ggcorrplot(cell_type_corr, method = "circle")


