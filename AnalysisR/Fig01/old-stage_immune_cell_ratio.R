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

# ROI information
roi_normal_lst <- list()
roi_aah_lst <- list()
roi_ais_lst <- list()
roi_mia_lst <- list()
roi_adc_lst <- list()
for (ind in 1:nrow(roi_df)) {
    if (roi_df$ROI_Location[ind] == "Normal" | roi_df$ROI_Location[ind] == "DistantNormal")
        roi_normal_lst <- append(roi_normal_lst, roi_df$ROI_ID[ind])
    else if (roi_df$ROI_Location[ind] == "Tumor") {
        if (roi_df$ROI_Diag[ind] == "AAH")
            roi_aah_lst <- append(roi_aah_lst, roi_df$ROI_ID[ind])
        else if (roi_df$ROI_Diag[ind] == "AIS")
            roi_ais_lst <- append(roi_ais_lst, roi_df$ROI_ID[ind])
        else if (roi_df$ROI_Diag[ind] == "MIA")
            roi_mia_lst <- append(roi_mia_lst, roi_df$ROI_ID[ind])
        else if (roi_df$ROI_Diag[ind] == "ADC")
            roi_adc_lst <- append(roi_adc_lst, roi_df$ROI_ID[ind])
        else
            print(roi_df$ROI_Diag[ind])
    }
}

## obtain cell-id & cell-type
cell_id_lst <- rownames(colData(spe))
celtype_lst <- spe$celltype

immune_lst <- c("CD4-T-Cell", "CD8-T-Cell", "T-Reg-Cell", "B-Cell", "Macrophage", 
                "Monocyte", "Dendritic-Cell", "Neutrophil", "MDSC", "NK-Cell")

# Normal ROI
normal_ratio_lst <- list()
for (ir in 1:length(roi_normal_lst)) {
    cur_roi = roi_normal_lst[[ir]]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- celtype_lst[cell_indices]
    ttl_num = length(cell_indices)
    immune_num = 0
    for (ii in 1:ttl_num) 
        if (roi_celltypes[ii] %in% immune_lst)
            immune_num <- immune_num + 1
    normal_ratio_lst <- append(normal_ratio_lst, immune_num * 1.0 / ttl_num)
}

# AAH ROI
aah_ratio_lst <- list()
for (ir in 1:length(roi_aah_lst)) {
    cur_roi = roi_aah_lst[[ir]]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- celtype_lst[cell_indices]
    ttl_num = length(cell_indices)
    immune_num = 0
    for (ii in 1:ttl_num) 
        if (roi_celltypes[ii] %in% immune_lst)
            immune_num <- immune_num + 1
    aah_ratio_lst <- append(aah_ratio_lst, immune_num * 1.0 / ttl_num)
}


# AIS ROI
ais_ratio_lst <- list()
for (ir in 1:length(roi_ais_lst)) {
    cur_roi = roi_ais_lst[[ir]]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- celtype_lst[cell_indices]
    ttl_num = length(cell_indices)
    immune_num = 0
    for (ii in 1:ttl_num) 
        if (roi_celltypes[ii] %in% immune_lst)
            immune_num <- immune_num + 1
    ais_ratio_lst <- append(ais_ratio_lst, immune_num * 1.0 / ttl_num)
}

# MIA ROI
mia_ratio_lst <- list()
for (ir in 1:length(roi_mia_lst)) {
    cur_roi = roi_mia_lst[[ir]]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- celtype_lst[cell_indices]
    ttl_num = length(cell_indices)
    immune_num = 0
    for (ii in 1:ttl_num) 
        if (roi_celltypes[ii] %in% immune_lst)
            immune_num <- immune_num + 1
    mia_ratio_lst <- append(mia_ratio_lst, immune_num * 1.0 / ttl_num)
}

# ADC ROI
adc_ratio_lst <- list()
for (ir in 1:length(roi_adc_lst)) {
    cur_roi = roi_adc_lst[[ir]]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- celtype_lst[cell_indices]
    ttl_num = length(cell_indices)
    immune_num = 0
    for (ii in 1:ttl_num) 
        if (roi_celltypes[ii] %in% immune_lst)
            immune_num <- immune_num + 1
    adc_ratio_lst <- append(adc_ratio_lst, immune_num * 1.0 / ttl_num)
}

roi_lst_all <- c(roi_normal_lst, roi_aah_lst, roi_ais_lst, roi_mia_lst, roi_adc_lst)
roi_stage_all <- c(rep("Normal", length(roi_normal_lst)), rep("AAH", length(roi_aah_lst)),
                   rep("AIS", length(roi_ais_lst)), rep("MIA", length(roi_mia_lst)), 
                   rep("ADC", length(roi_adc_lst)))
roi_ratio_all <- c(normal_ratio_lst, aah_ratio_lst, ais_ratio_lst, mia_ratio_lst, adc_ratio_lst)
stage_roi_ratio_df <- data.frame(ID = unlist(roi_lst_all), 
                                 Stage = roi_stage_all,
                                 Ratio = unlist(roi_ratio_all))

MinMeanSEMMax <- function(x) {
    v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
}

stage_order <- c("Normal", "AAH", "AIS", "MIA", "ADC")
ggplot(stage_roi_ratio_df, aes(x = factor(Stage, level=stage_order), y=Ratio)) + 
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="black") + 
    geom_beeswarm(size=1.0) + ggtitle("Beeswarm")
