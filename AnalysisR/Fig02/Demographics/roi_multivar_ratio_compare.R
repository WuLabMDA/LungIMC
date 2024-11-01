library(imcRtools)
library(readxl)
library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(gridExtra)
library(dplyr)

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


roi_lst <- roi_df$ROI_ID
roi_vec <- c()
for (ind in 1:nrow(roi_df)) {
    if (roi_df$ROI_Location[ind] == "Tumor" ) 
        roi_vec <- append(roi_vec, roi_df$ROI_ID[ind])
}
roi_slide_lst <- str_extract_all(roi_vec, ".+(?=-ROI)", simplify = TRUE)

slide_info_name <- "Lesion_Info"
slide_info_path <- file.path(metadata_dir, paste0(slide_info_name, ".xlsx"))
slide_df <- read_excel(slide_info_path)
# pathological stages
# AAH/AIS/MIA/ADC
path_stage <- "ADC"
stage_slide_df <- filter(slide_df, Slide_Diag == path_stage)
stage_slide_lst <- stage_slide_df$Slide_ID

## obtain lymphoid ratio
cell_id_lst <- rownames(colData(spe))
celtype_lst <- spe$celltype

# load patient information
patient_info_name <- "Patient_Info"
patient_info_path <- file.path(metadata_dir, paste0(patient_info_name, ".xlsx"))
patient_df <- read_excel(patient_info_path)

# MinMeanSEMMax
MinMeanSEMMax <- function(x) {
    v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
}

# set the cell-type
cell_type <- "Proliferating-Cell"

roi_lst <- c()
ratio_lst <- c()
gender_lst <- c()
race_lst <- c()
median_age <- 71.3
age_lst <- c()
smoke_lst <- c()

for (ir in 1:length(roi_vec)) {
    cur_roi <- roi_vec[ir]
    cur_slide <- roi_slide_lst[ir]
    if (!(cur_slide %in% stage_slide_lst))
        next
    roi_lst <- append(roi_lst, cur_roi)
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- celtype_lst[cell_indices]
    ttl_num <- length(roi_celltypes)
    cell_num <- sum(roi_celltypes == cell_type)
    cell_ratio <- cell_num * 1.0 / ttl_num
    ratio_lst <- append(ratio_lst, cell_ratio)
    underline_pos <- unlist(gregexpr("-", cur_roi))
    pat_id <- substr(cur_roi, 1, underline_pos[length(underline_pos)-1]-1)
    df_index <- which(patient_df$PatientID == pat_id)
    gender_lst <- append(gender_lst, patient_df$Gender[df_index])
    race_lst <- append(race_lst, patient_df$Race[df_index])
    
    if (patient_df$Age[df_index] <= median_age) 
        age_lst <- append(age_lst, "<=71")
    else
        age_lst <- append(age_lst, ">71")
    smoke_lst <- append(smoke_lst, patient_df$SmokeStatus[df_index])
}

## Construct data frame
cell_ratio_df <- data.frame(ROI=roi_lst, Ratio=ratio_lst, Gender=gender_lst, Race=race_lst,
                            Age=age_lst, Smoke=smoke_lst)
## Plotting
p_gender <- ggplot(cell_ratio_df, aes(x = factor(Gender, level=c("M", "F")), y=Ratio)) + 
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="black") + 
    geom_beeswarm(cex = 2.5, corral = "random", corral.width = 0.5)

p_race <- ggplot(cell_ratio_df, aes(x = factor(Race, level=c("Asian", "White")), y=Ratio)) + 
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="black") + 
    geom_beeswarm(cex = 2.5, corral = "random", corral.width = 0.5)

p_age <- ggplot(cell_ratio_df, aes(x = factor(Age, level=c("<=71", ">71")), y=Ratio)) + 
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="black") + 
    geom_beeswarm(cex = 2.5, corral = "random", corral.width = 0.5)

p_smoke <- ggplot(cell_ratio_df, aes(x = factor(Smoke, level=c("Non-Smoker", "Smoker")), y=Ratio)) + 
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="black") + 
    geom_beeswarm(cex = 2.5, corral = "random", corral.width = 0.5)

grid.arrange(p_gender, p_race, p_age, p_smoke, ncol=4)
