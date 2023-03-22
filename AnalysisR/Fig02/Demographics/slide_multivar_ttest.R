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
    roi_location <- roi_df$ROI_Location[ind]
    if (roi_location %in% c("Normal", "Tumor")) 
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


all_cell_lst <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell",
                  "CD8-T-Cell", "CD4-T-Cell",  "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte",
                  "MDSC", "Fibroblast", "Undefined")
age_pvals <- c()
age_cmps <- c()
gender_pvals <- c()
gender_cmps <- c()
race_pvals <- c()
race_cmps <- c()
# recur_pvals <- c()
smoke_pvals <- c()
smoke_cmps <- c()
normal_aah_pvals <- c()
normal_aah_cmps <- c()
normal_ais_pvals <- c()
normal_ais_cmps <- c()
normal_mia_pvals <- c()
normal_mia_cmps <- c()
normal_adc_pvals <- c()
normal_adc_cmps <- c()
aah_ais_pvals <- c()
aah_ais_cmps <- c()
aah_mia_pvals <- c()
aah_mia_cmps <- c()
aah_adc_pvals <- c()
aah_adc_cmps <- c()
ais_mia_pvals <- c()
ais_mia_cmps <- c()
ais_adc_pvals <- c()
ais_adc_cmps <- c()
mia_adc_pvals <- c()
mia_adc_cmps <- c()

for (cell_type in all_cell_lst) {
    # collect information
    ratio_lst <- c()
    gender_lst <- c()
    race_lst <- c()
    # recurrent_lst <- c()
    median_age <- 71.3
    age_lst <- c()
    smoke_lst <- c()
    stage_lst <- c()
    
    # iterate by slide
    stage_slide_lst <- slide_info_df$Slide_ID
    for (ir in 1:length(stage_slide_lst)) {
        cur_slide <- stage_slide_lst[ir]
        if (cur_slide == "2571-1D")
            next
        # if (cur_slide == "H12-0330-2")
        #     next    
        slide_cell_indices <- which(startsWith(cell_id_lst, cur_slide))
        slide_celltypes <- celltype_lst[slide_cell_indices]
        slide_cell_ratio <- sum(slide_celltypes == cell_type) / length(slide_celltypes)
        ratio_lst <- append(ratio_lst, slide_cell_ratio)
        underline_pos <- unlist(gregexpr("-", cur_slide))
        pat_id <- substr(cur_slide, 1, underline_pos[length(underline_pos)]-1)
        patient_index <- which(patient_df$PatientID == pat_id)
        gender_lst <- append(gender_lst, patient_df$Gender[patient_index])
        race_lst <- append(race_lst, patient_df$Race[patient_index])
        # if (patient_df$RecurrentStatus[patient_index] == "RECURRENT")
        #     recurrent_lst <- append(recurrent_lst, "Recur")
        # else
        #     recurrent_lst <- append(recurrent_lst, "Non-Recur")
        if (patient_df$Age[patient_index] <= median_age) 
            age_lst <- append(age_lst, "<=71")
        else
            age_lst <- append(age_lst, ">71")
        smoke_lst <- append(smoke_lst, patient_df$SmokeStatus[patient_index])
        slide_index <- which(slide_info_df$Slide_ID == cur_slide)
        # if (length(slide_index) == 2)
        #     print(cur_slide)
        stage_lst <- append(stage_lst, slide_info_df$Slide_Diag[slide_index])
    }
    
    cell_ratio_df <- data.frame(Ratio=ratio_lst, Gender=gender_lst, Race=race_lst, Age=age_lst, Smoke=smoke_lst, Stage=stage_lst)
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio)
    
    # Welch Two Sample t-test p-value
    age_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="<=71"], long_ratio_df$Ratio[long_ratio_df$value==">71"])
    age_pvals <- append(age_pvals, age_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value==">71"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="<=71"]))
        age_cmps <- append(age_cmps, TRUE)
    else
        age_cmps <- append(age_cmps, FALSE)
    gender_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="F"], long_ratio_df$Ratio[long_ratio_df$value=="M"])
    gender_pvals <- append(gender_pvals, gender_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="F"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="M"]))
        gender_cmps <- append(gender_cmps, TRUE)
    else
        gender_cmps <- append(gender_cmps, FALSE)
    race_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="White"], long_ratio_df$Ratio[long_ratio_df$value=="Asian"])
    race_pvals <- append(race_pvals, race_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="White"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Asian"]))
        race_cmps <- append(race_cmps, TRUE)
    else
        race_cmps <- append(race_cmps, FALSE)    
    # recurrent_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Non-Recur"], long_ratio_df$Ratio[long_ratio_df$value=="Recur"])
    # recur_pvals <- append(recur_pvals, recurrent_ttest$p.value)
    smoke_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Non-Smoker"], long_ratio_df$Ratio[long_ratio_df$value=="Smoker"])
    smoke_pvals <- append(smoke_pvals, smoke_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Smoker"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Non-Smoker"]))
        smoke_cmps <- append(smoke_cmps, TRUE)
    else
        smoke_cmps <- append(smoke_cmps, FALSE)  
    normal_aah_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Normal"], long_ratio_df$Ratio[long_ratio_df$value=="AAH"])
    normal_aah_pvals <- append(normal_aah_pvals, normal_aah_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="AAH"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Normal"]))
        normal_aah_cmps <- append(normal_aah_cmps, TRUE)
    else
        normal_aah_cmps <- append(normal_aah_cmps, FALSE)     
    normal_ais_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Normal"], long_ratio_df$Ratio[long_ratio_df$value=="AIS"])
    normal_ais_pvals <- append(normal_ais_pvals, normal_ais_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="AIS"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Normal"]))
        normal_ais_cmps <- append(normal_ais_cmps, TRUE)
    else
        normal_ais_cmps <- append(normal_ais_cmps, FALSE)      
    normal_mia_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Normal"], long_ratio_df$Ratio[long_ratio_df$value=="MIA"])
    normal_mia_pvals <- append(normal_mia_pvals, normal_mia_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="MIA"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Normal"]))
        normal_mia_cmps <- append(normal_mia_cmps, TRUE)
    else
        normal_mia_cmps <- append(normal_mia_cmps, FALSE)       
    normal_adc_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Normal"], long_ratio_df$Ratio[long_ratio_df$value=="ADC"])
    normal_adc_pvals <- append(normal_adc_pvals, normal_adc_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="ADC"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Normal"]))
        normal_adc_cmps <- append(normal_adc_cmps, TRUE)
    else
        normal_adc_cmps <- append(normal_adc_cmps, FALSE)  
    aah_ais_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AIS"], long_ratio_df$Ratio[long_ratio_df$value=="AAH"])
    aah_ais_pvals <- append(aah_ais_pvals, aah_ais_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="AIS"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AAH"]))
        aah_ais_cmps <- append(aah_ais_cmps, TRUE)
    else
        aah_ais_cmps <- append(aah_ais_cmps, FALSE) 
    aah_mia_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AAH"], long_ratio_df$Ratio[long_ratio_df$value=="MIA"])
    aah_mia_pvals <- append(aah_mia_pvals, aah_mia_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="MIA"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AAH"]))
        aah_mia_cmps <- append(aah_mia_cmps, TRUE)
    else
        aah_mia_cmps <- append(aah_mia_cmps, FALSE)     
    aah_adc_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AAH"], long_ratio_df$Ratio[long_ratio_df$value=="ADC"])
    aah_adc_pvals <- append(aah_adc_pvals, aah_adc_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="ADC"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AAH"]))
        aah_adc_cmps <- append(aah_adc_cmps, TRUE)
    else
        aah_adc_cmps <- append(aah_adc_cmps, FALSE)      
    ais_mia_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AIS"], long_ratio_df$Ratio[long_ratio_df$value=="MIA"])
    ais_mia_pvals <- append(ais_mia_pvals, ais_mia_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="MIA"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AIS"]))
        ais_mia_cmps <- append(ais_mia_cmps, TRUE)
    else
        ais_mia_cmps <- append(ais_mia_cmps, FALSE)       
    ais_adc_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AIS"], long_ratio_df$Ratio[long_ratio_df$value=="ADC"])
    ais_adc_pvals <- append(ais_adc_pvals, ais_adc_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="ADC"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AIS"]))
        ais_adc_cmps <- append(ais_adc_cmps, TRUE)
    else
        ais_adc_cmps <- append(ais_adc_cmps, FALSE)       
    mia_adc_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="MIA"], long_ratio_df$Ratio[long_ratio_df$value=="ADC"])
    mia_adc_pvals <- append(mia_adc_pvals, mia_adc_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="ADC"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="MIA"]))
        mia_adc_cmps <- append(mia_adc_cmps, TRUE)
    else
        mia_adc_cmps <- append(mia_adc_cmps, FALSE)     
    
}

p_val_df <- data.frame(Age=age_pvals, Gender=gender_pvals, Race=race_pvals, Smoke=smoke_pvals,
                       Normal_AAH=normal_aah_pvals, Normal_AIS=normal_ais_pvals, Normal_MIA=normal_mia_pvals, Normal_ADC=normal_adc_pvals,
                       AAH_AIS=aah_ais_pvals, AAH_MIA=aah_mia_pvals, AAH_ADC=aah_adc_pvals, AIS_MIA=ais_mia_pvals,
                       AIS_ADC=ais_adc_pvals, MIA_ADC=mia_adc_pvals, row.names = all_cell_lst)
p_cmp_df <- data.frame(Age=age_cmps, Gender=gender_cmps, Race=race_cmps, Smoke=smoke_cmps,
                       Normal_AAH=normal_aah_cmps, Normal_AIS=normal_ais_cmps, Normal_MIA=normal_mia_cmps, Normal_ADC=normal_adc_cmps,
                       AAH_AIS=aah_ais_cmps, AAH_MIA=aah_mia_cmps, AAH_ADC=aah_adc_cmps, AIS_MIA=ais_mia_cmps,
                       AIS_ADC=ais_adc_cmps, MIA_ADC=mia_adc_cmps, row.names = all_cell_lst)


# var_order <- c("Age(>71)", "Gender(F)", "Race(White)", "Non-Smoker>Smoker", "AAH>Normal", "AIS>Normal", "MIA>Normal",
#                "ADC>Normal", "AIS>AAH", "MIA>AAH", "ADC>AAH", "MIA>AIS", "ADC>AIS", "ADC>MIA")
var_order <- c("Age", "Gender", "Race", "Smoke", "Normal_AAH", "Normal_AIS", "Normal_MIA",
               "Normal_ADC", "AAH_AIS", "AAH_MIA", "AAH_ADC", "AIS_MIA", "AIS_ADC", "MIA_ADC")

group_df <- p_val_df %>% rownames_to_column(var = "CellType") %>% gather(Vars, Pvals, -CellType) 
group_df$AdjustPs <- p.adjust(group_df$Pvals, method = "fdr")
# adjust_df <- group_df %>% mutate(gGroup = case_when(AdjustPs > 0.05 ~ 'NS', AdjustPs > 0.01 ~ '*', AdjustPs > 0.001 ~ '**', .default = "***"))
adjust_df <- group_df %>% mutate(gGroup = case_when(AdjustPs > 0.05 ~ 'NS', AdjustPs > 0.01 ~ '*', .default = "**"))

cmp_df <- p_cmp_df %>% rownames_to_column(var = "CellType") %>% gather(Vars, Cmps, -CellType) 
adjust_df$Cmps <- as.factor(cmp_df$Cmps) 

adjust_df %>% ggplot() +
    geom_point(aes(x = factor(Vars, level=var_order), y = factor(CellType, level=rev(all_cell_lst)), 
                   size = p.to.Z(Pvals), shape = factor(gGroup, level=c("NS", "*", "**")),
                   color = Cmps)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 



