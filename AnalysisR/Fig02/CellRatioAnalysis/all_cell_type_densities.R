library(imcRtools)
library(readxl)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

metadata_dir <- file.path(data_root_dir, "Metadata")
# load ROI information
roi_info_name <- "ROI_Info_Aggregation.csv"
roi_info_path <- file.path(metadata_dir, roi_info_name)
roi_df <- read_csv(roi_info_path)

# ROI list
roi_vec <- roi_df$ROI_ID
roi_area <- roi_df$Area


## obtain lymphoid ratio
cell_id_lst <- rownames(colData(spe))
celltype_lst <- spe$celltype


all_cell_lst <- c("Epithelial-Cell", "Endothelial-Cell", "Fibroblast", "CD4-T-Cell", "CD8-T-Cell", 
                  "T-Reg-Cell", "B-Cell", "Macrophage", "Monocyte", "Dendritic-Cell", 
                  "Neutrophil", "MDSC", "NK-Cell", "Proliferating-Cell", "Undefined")
all_type_lst <- c()
all_ratio_lst <- c()

for (ir in 1:length(roi_vec)) {
    cur_roi = roi_vec[ir]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- celltype_lst[cell_indices]
    for (cell_type in all_cell_lst) {
        all_type_lst <- append(all_type_lst, cell_type)
        all_ratio_lst <- append(all_ratio_lst, sum(roi_celltypes == cell_type) * 1000000 / roi_area[ir])
    }
}

celltype_mean_ratio_lst <- list()
for (cell_type in all_cell_lst) {
    cell_indices <- which(all_type_lst == cell_type)
    cell_ratio <- mean(all_ratio_lst[cell_indices])
    celltype_mean_ratio_lst <- append(celltype_mean_ratio_lst, cell_ratio)
}
names(celltype_mean_ratio_lst) <- all_cell_lst
celltype_mean_ratio_lst <- celltype_mean_ratio_lst[order(unlist(celltype_mean_ratio_lst), decreasing=TRUE)]

# Construct data frame
cell_type_density_df <- data.frame(Type=all_type_lst, Density=all_ratio_lst)
for (cell_type in all_cell_lst) {
    print(cell_type)
    ct_df <- filter(cell_type_density_df, Type == cell_type)
    ct_ratios <- ct_df$Density
    ct_mean <- mean(ct_ratios)
    ct_sd <- sd(ct_ratios)
    print(paste0("Mean: ", ct_mean, "SD: ", ct_sd))
}

MinMeanSEMMax <- function(x) {
    v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
}
ggplot(cell_type_density_df, aes(x = factor(Type, level=names(celltype_mean_ratio_lst)), y=Density)) + 
    theme(axis.text.x=element_text(angle=90, hjust=1)) + 
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="black") + 
    geom_beeswarm(cex = 2.5, corral = "random", corral.width = 0.4)
