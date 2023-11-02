library(imcRtools)
library(readxl)
library(ggplot2)
library(ggbeeswarm)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_name <- "ROI_Info"
roi_info_path <- file.path(metadata_dir, paste0(roi_info_name, ".xlsx"))
roi_df <- read_excel(roi_info_path)

# ROI list
roi_lst <- roi_df$ROI_ID

# ROI stage information
roi_vec <- c()
stage_vec <- c()
for (ind in 1:nrow(roi_df)) {
    if (roi_df$ROI_Location[ind] == "Normal" | roi_df$ROI_Location[ind] == "DistantNormal") {
        roi_vec <- append(roi_vec, roi_df$ROI_ID[ind])
        stage_vec <- append(stage_vec, "Normal")        
    }
    else if (roi_df$ROI_Location[ind] == "Tumor") {
        roi_vec <- append(roi_vec, roi_df$ROI_ID[ind])
        stage_vec <- append(stage_vec, roi_df$ROI_Diag[ind])
    }
}

## obtain myeloid ratio
cell_id_lst <- rownames(colData(spe))
celtype_lst <- spe$celltype

# all_cell_lst <- c("Epithelial-Cell", "Endothelial-Cell", "Fibroblast", "CD4-T-Cell", "CD8-T-Cell", 
#                 "T-Reg-Cell", "B-Cell", "Macrophage", "Monocyte", "Dendritic-Cell", 
#                 "Neutrophil", "MDSC", "NK-Cell", "Proliferating-Cell", "Undefined")

interested_celltype <- "Proliferating-Cell"
cell_ratio_vec <- c()
for (ir in 1:length(roi_vec)) {
    cur_roi = roi_vec[ir]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- celtype_lst[cell_indices]
    ttl_num = length(roi_celltypes)
    cell_num <- sum(roi_celltypes == interested_celltype)
    cell_ratio_vec <- append(cell_ratio_vec, cell_num * 1.0 / ttl_num)
}

# Construct data frame
stage_cell_ratio_df <- data.frame(ROI = roi_vec, Stage = stage_vec, Ratio = cell_ratio_vec)

MinMeanSEMMax <- function(x) {
    v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
}

stage_order <- c("Normal", "AAH", "AIS", "MIA", "ADC")
ggplot(stage_cell_ratio_df, aes(x = factor(Stage, level=stage_order), y=Ratio)) + 
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="black") + 
    geom_beeswarm(cex = 2.5, corral = "random", corral.width = 0.4)
