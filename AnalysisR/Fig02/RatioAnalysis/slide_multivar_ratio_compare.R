library(imcRtools)
library(readxl)
library(stringr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(gridExtra)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

metadata_dir <- file.path(data_root_dir, "Metadata")
# load ROI information
roi_info_name <- "ROI_Info"
roi_info_path <- file.path(metadata_dir, paste0(roi_info_name, ".xlsx"))
roi_df <- read_excel(roi_info_path)
# load lesion information
slide_info_name <- "Lesion_Info"
slide_info_path <- file.path(metadata_dir, paste0(slide_info_name, ".xlsx"))
slide_info_df <- read_excel(slide_info_path)
# load patient information
patient_info_name <- "Patient_Info"
patient_info_path <- file.path(metadata_dir, paste0(patient_info_name, ".xlsx"))
patient_df <- read_excel(patient_info_path)

## filtering rois
roi_lst <- roi_df$ROI_ID
roi_vec <- c()
for (ind in 1:nrow(roi_df)) {
    if (roi_df$ROI_Location[ind] == "Tumor" ) 
        roi_vec <- append(roi_vec, roi_df$ROI_ID[ind])
}
## obtain cells inside tumor lesions
lesion_indices = c()
cell_id_lst <- rownames(colData(spe))
for (cur_roi in roi_vec) 
    lesion_indices <- append(lesion_indices, which(startsWith(cell_id_lst, cur_roi)))
cell_id_lst <- cell_id_lst[lesion_indices]
celltype_lst <- spe$celltype
celltype_lst <- celltype_lst[lesion_indices]


# set pathological stages (AAH/AIS/MIA/ADC)
path_stage <- "ADC"
cell_type <- "Proliferating-Cell"

# filtering stages
stage_slide_df <- filter(slide_info_df, Slide_Diag==path_stage)
stage_slide_lst <- stage_slide_df$Slide_ID


ratio_lst <- c()
gender_lst <- c()
race_lst <- c()
recurrent_lst <- c()
median_age <- 71.3
age_lst <- c()
smoke_lst <- c()

for (ir in 1:length(stage_slide_lst)) {
    cur_slide <- stage_slide_lst[ir]
    if (cur_slide == "H12-0330-2")
        next
    slide_cell_indices <- which(startsWith(cell_id_lst, cur_slide))
    slide_celltypes <- celltype_lst[slide_cell_indices]
    slide_cell_ratio <- sum(slide_celltypes == cell_type) / length(slide_celltypes)
    ratio_lst <- append(ratio_lst, slide_cell_ratio)
    underline_pos <- unlist(gregexpr("-", cur_slide))
    pat_id <- substr(cur_slide, 1, underline_pos[length(underline_pos)]-1)
    df_index <- which(patient_df$PatientID == pat_id)
    gender_lst <- append(gender_lst, patient_df$Gender[df_index])
    race_lst <- append(race_lst, patient_df$Race[df_index])
    if (patient_df$RecurrentStatus[df_index] == "RECURRENT")
        recurrent_lst <- append(recurrent_lst, "Recur")
    else
        recurrent_lst <- append(recurrent_lst, "Non-Recur")
    if (patient_df$Age[df_index] <= median_age) 
        age_lst <- append(age_lst, "<=71")
    else
        age_lst <- append(age_lst, ">71")
    smoke_lst <- append(smoke_lst, patient_df$SmokeStatus[df_index])
}

# MinMeanSEMMax
MinMeanSEMMax <- function(x) {
    v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
}


cell_ratio_df <- data.frame(Ratio=ratio_lst, Gender=gender_lst, Race=race_lst,
                            Recurrent=recurrent_lst, Age=age_lst, Smoke=smoke_lst)
long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio)
ggplot(long_ratio_df, aes(x = value, y = Ratio) ) +
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="black") + 
    geom_beeswarm(cex = 2.5, corral = "random", corral.width = 0.5) +
    facet_wrap(~variable, scales = "free_x", ncol = 5)
