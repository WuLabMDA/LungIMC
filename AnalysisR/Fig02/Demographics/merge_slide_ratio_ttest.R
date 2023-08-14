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
roi_info_name <- "ROI_Info_Aggregation.csv"
roi_info_path <- file.path(metadata_dir, roi_info_name)
roi_df <- read_csv(roi_info_path)
# load lesion information
slide_info_name <- "Lesion_Info_AggregationSize.csv"
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
    if (roi_location %in% c("Normal", "Tumor")) 
        interested_roi_vec <- append(interested_roi_vec, roi_df$ROI_ID[ind])
}

## obtain cells inside tumor lesions
lesion_indices = c()
cell_id_lst <- rownames(colData(spe))
for (cur_roi in interested_roi_vec) 
    lesion_indices <- append(lesion_indices, which(startsWith(cell_id_lst, cur_roi)))
cell_id_lst <- cell_id_lst[lesion_indices]
celltype_lst <- spe$celltype
celltype_lst <- celltype_lst[lesion_indices]

# list all cell types
all_cell_lst <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell",
                  "CD8-T-Cell", "CD4-T-Cell",  "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte",
                  "MDSC", "Fibroblast", "Undefined")

age_pvals <- c()
age_cmps <- c()
gender_pvals <- c()
gender_cmps <- c()
race_pvals <- c()
race_cmps <- c()
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
aah_tmb_pvals <- c()
aah_tmb_cmps <- c()
ais_tmb_pvals <- c()
ais_tmb_cmps <- c()
mia_tmb_pvals <- c()
mia_tmb_cmps <- c()
adc_tmb_pvals <- c()
adc_tmb_cmps <- c()
mia_egfr_kras_pvals <- c()
mia_egfr_kras_cmps <- c()
mia_egfr_other_pvals <- c()
mia_egfr_other_cmps <- c()
mia_kras_other_pvals <- c()
mia_kras_other_cmps <- c()
adc_egfr_kras_pvals <- c()
adc_egfr_kras_cmps <- c()
adc_egfr_other_pvals <- c()
adc_egfr_other_cmps <- c()
adc_kras_other_pvals <- c()
adc_kras_other_cmps <- c()
aah_size_pvals <- c()
aah_size_cmps <- c()
ais_size_pvals <- c()
ais_size_cmps <- c()
mia_size_pvals <- c()
mia_size_cmps <- c()
adc_size_pvals <- c()
adc_size_cmps <- c()


for (cell_type in all_cell_lst) {
    # collect information
    ratio_lst <- c()
    gender_lst <- c()
    race_lst <- c()
    median_age <- 71.3
    age_lst <- c()
    smoke_lst <- c()
    stage_lst <- c()
    
    aah_tmb_ratio_lst <- c()
    aah_tmb_status_lst <- c()
    ais_tmb_ratio_lst <- c()
    ais_tmb_status_lst <- c()
    mia_tmb_ratio_lst <- c()
    mia_tmb_status_lst <- c()
    adc_tmb_ratio_lst <- c()
    adc_tmb_status_lst <- c()
    
    mia_egfr_kras_ratio_lst <- c()
    mia_egfr_kras_status_lst <- c()
    adc_egfr_kras_ratio_lst <- c()
    adc_egfr_kras_status_lst <- c()
    
    aah_size_ratio_lst <- c()
    aah_size_status_lst <- c()
    ais_size_ratio_lst <- c()
    ais_size_status_lst <- c()
    mia_size_ratio_lst <- c()
    mia_size_status_lst <- c()
    adc_size_ratio_lst <- c()
    adc_size_status_lst <- c()    
    
    # iterate by slide
    stage_slide_lst <- slide_info_df$Slide_ID
    for (ir in 1:length(stage_slide_lst)) {
        # iterate by slide
        cur_slide <- stage_slide_lst[ir]
        if (cur_slide == "2571-1D")
            next
        slide_cell_indices <- which(startsWith(cell_id_lst, cur_slide))
        slide_celltypes <- celltype_lst[slide_cell_indices]
        slide_cell_ratio <- sum(slide_celltypes == cell_type) / length(slide_celltypes)
        
        # gather information
        ratio_lst <- append(ratio_lst, slide_cell_ratio)
        underline_pos <- unlist(gregexpr("-", cur_slide))
        pat_id <- substr(cur_slide, 1, underline_pos[length(underline_pos)]-1)
        patient_index <- which(patient_df$PatientID == pat_id)
        gender_lst <- append(gender_lst, patient_df$Gender[patient_index])
        race_lst <- append(race_lst, patient_df$Race[patient_index])
        if (patient_df$Age[patient_index] <= median_age) 
            age_lst <- append(age_lst, "<=71")
        else
            age_lst <- append(age_lst, ">71")
        smoke_lst <- append(smoke_lst, patient_df$SmokeStatus[patient_index])
        slide_index <- which(slide_info_df$Slide_ID == cur_slide)
        stage_lst <- append(stage_lst, slide_info_df$Slide_Diag[slide_index])    
        
        cur_diag <- slide_info_df$Slide_Diag[ir]
        cur_tmb <- slide_info_df$TMB[ir]
        cur_egfr_kras <- slide_info_df$EGFR_KRAS[ir]
        cur_size <- slide_info_df$TumorSize[ir]
        if (!is.na(cur_tmb)) {
            if (cur_diag == "AAH") {
                aah_tmb_ratio_lst <- append(aah_tmb_ratio_lst, slide_cell_ratio)
                aah_tmb_status_lst <- append(aah_tmb_status_lst, cur_tmb)
            }
            if (cur_diag == "AIS") {
                ais_tmb_ratio_lst <- append(ais_tmb_ratio_lst, slide_cell_ratio)
                ais_tmb_status_lst <- append(ais_tmb_status_lst, cur_tmb)
            }
            if (cur_diag == "MIA") {
                mia_tmb_ratio_lst <- append(mia_tmb_ratio_lst, slide_cell_ratio)
                mia_tmb_status_lst <- append(mia_tmb_status_lst, cur_tmb)
            }
            if (cur_diag == "ADC") {
                adc_tmb_ratio_lst <- append(adc_tmb_ratio_lst, slide_cell_ratio)
                adc_tmb_status_lst <- append(adc_tmb_status_lst, cur_tmb)
            }            
        }
        if (!is.na(cur_egfr_kras)) {
            if (cur_diag == "MIA") {
                mia_egfr_kras_ratio_lst <- append(mia_egfr_kras_ratio_lst, slide_cell_ratio)
                mia_egfr_kras_status_lst <- append(mia_egfr_kras_status_lst, cur_egfr_kras)
            }
            if (cur_diag == "ADC") {
                adc_egfr_kras_ratio_lst <- append(adc_egfr_kras_ratio_lst, slide_cell_ratio)
                adc_egfr_kras_status_lst <- append(adc_egfr_kras_status_lst, cur_egfr_kras)
            }              
        }
        if (!is.na(cur_size)) {
            if (cur_diag == "AAH") {
                aah_size_ratio_lst <- append(aah_size_ratio_lst, slide_cell_ratio)
                aah_size_status_lst <- append(aah_size_status_lst, cur_size)
            }
            if (cur_diag == "AIS") {
                ais_size_ratio_lst <- append(ais_size_ratio_lst, slide_cell_ratio)
                ais_size_status_lst <- append(ais_size_status_lst, cur_size)
            }
            if (cur_diag == "MIA") {
                mia_size_ratio_lst <- append(mia_size_ratio_lst, slide_cell_ratio)
                mia_size_status_lst <- append(mia_size_status_lst, cur_size)
            }
            if (cur_diag == "ADC") {
                adc_size_ratio_lst <- append(adc_size_ratio_lst, slide_cell_ratio)
                adc_size_status_lst <- append(adc_size_status_lst, cur_size)
            }            
        }        
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
    
    
    ## AAH TMB
    cell_ratio_df <- data.frame(Ratio=aah_tmb_ratio_lst, TMB=aah_tmb_status_lst)    
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 
    aah_tmb_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Low"], long_ratio_df$Ratio[long_ratio_df$value=="High"])
    aah_tmb_pvals <- append(aah_tmb_pvals, aah_tmb_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="High"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Low"]))
        aah_tmb_cmps <- append(aah_tmb_cmps, TRUE)
    else
        aah_tmb_cmps <- append(aah_tmb_cmps, FALSE)      
    ## AIS TMB
    cell_ratio_df <- data.frame(Ratio=ais_tmb_ratio_lst, TMB=ais_tmb_status_lst)    
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 
    ais_tmb_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Low"], long_ratio_df$Ratio[long_ratio_df$value=="High"])
    ais_tmb_pvals <- append(ais_tmb_pvals, ais_tmb_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="High"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Low"]))
        ais_tmb_cmps <- append(ais_tmb_cmps, TRUE)
    else
        ais_tmb_cmps <- append(ais_tmb_cmps, FALSE)  
    ## MIA TMB
    cell_ratio_df <- data.frame(Ratio=mia_tmb_ratio_lst, TMB=mia_tmb_status_lst)    
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 
    mia_tmb_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Low"], long_ratio_df$Ratio[long_ratio_df$value=="High"])
    mia_tmb_pvals <- append(mia_tmb_pvals, mia_tmb_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="High"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Low"]))
        mia_tmb_cmps <- append(mia_tmb_cmps, TRUE)
    else
        mia_tmb_cmps <- append(mia_tmb_cmps, FALSE) 
    ## ADC TMB
    cell_ratio_df <- data.frame(Ratio=adc_tmb_ratio_lst, TMB=adc_tmb_status_lst)    
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 
    adc_tmb_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Low"], long_ratio_df$Ratio[long_ratio_df$value=="High"])
    adc_tmb_pvals <- append(adc_tmb_pvals, adc_tmb_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="High"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Low"]))
        adc_tmb_cmps <- append(adc_tmb_cmps, TRUE)
    else
        adc_tmb_cmps <- append(adc_tmb_cmps, FALSE)  
    
    ## MIA Mutation
    cell_ratio_df <- data.frame(Ratio=mia_egfr_kras_ratio_lst, EGFR_KRAS=mia_egfr_kras_status_lst)    
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 
    mia_egfr_kras_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="EGFR"], long_ratio_df$Ratio[long_ratio_df$value=="KRAS"])
    mia_egfr_kras_pvals <- append(mia_egfr_kras_pvals, mia_egfr_kras_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="EGFR"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="KRAS"]))
        mia_egfr_kras_cmps <- append(mia_egfr_kras_cmps, TRUE)
    else
        mia_egfr_kras_cmps <- append(mia_egfr_kras_cmps, FALSE)      
    mia_egfr_other_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="EGFR"], long_ratio_df$Ratio[long_ratio_df$value=="Other"])
    mia_egfr_other_pvals <- append(mia_egfr_other_pvals, mia_egfr_other_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="EGFR"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Other"]))
        mia_egfr_other_cmps <- append(mia_egfr_other_cmps, TRUE)
    else
        mia_egfr_other_cmps <- append(mia_egfr_other_cmps, FALSE)  
    mia_kras_other_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="KRAS"], long_ratio_df$Ratio[long_ratio_df$value=="Other"])
    mia_kras_other_pvals <- append(mia_kras_other_pvals, mia_kras_other_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="KRAS"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Other"]))
        mia_kras_other_cmps <- append(mia_kras_other_cmps, TRUE)
    else
        mia_kras_other_cmps <- append(mia_kras_other_cmps, FALSE)   
    
    ## MIA Mutation
    cell_ratio_df <- data.frame(Ratio=adc_egfr_kras_ratio_lst, EGFR_KRAS=adc_egfr_kras_status_lst)    
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 
    adc_egfr_kras_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="EGFR"], long_ratio_df$Ratio[long_ratio_df$value=="KRAS"])
    adc_egfr_kras_pvals <- append(adc_egfr_kras_pvals, adc_egfr_kras_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="EGFR"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="KRAS"]))
        adc_egfr_kras_cmps <- append(adc_egfr_kras_cmps, TRUE)
    else
        adc_egfr_kras_cmps <- append(adc_egfr_kras_cmps, FALSE)      
    adc_egfr_other_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="EGFR"], long_ratio_df$Ratio[long_ratio_df$value=="Other"])
    adc_egfr_other_pvals <- append(adc_egfr_other_pvals, adc_egfr_other_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="EGFR"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Other"]))
        adc_egfr_other_cmps <- append(adc_egfr_other_cmps, TRUE)
    else
        adc_egfr_other_cmps <- append(adc_egfr_other_cmps, FALSE)  
    adc_kras_other_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="KRAS"], long_ratio_df$Ratio[long_ratio_df$value=="Other"])
    adc_kras_other_pvals <- append(adc_kras_other_pvals, adc_kras_other_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="KRAS"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Other"]))
        adc_kras_other_cmps <- append(adc_kras_other_cmps, TRUE)
    else
        adc_kras_other_cmps <- append(adc_kras_other_cmps, FALSE)    
    
    ## AAH Size
    cell_ratio_df <- data.frame(Ratio=aah_size_ratio_lst, Size=aah_size_status_lst)    
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 
    aah_size_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Low"], long_ratio_df$Ratio[long_ratio_df$value=="High"])
    aah_size_pvals <- append(aah_size_pvals, aah_size_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="High"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Low"]))
        aah_size_cmps <- append(aah_size_cmps, TRUE)
    else
        aah_size_cmps <- append(aah_size_cmps, FALSE)      
    ## AIS Size
    cell_ratio_df <- data.frame(Ratio=ais_size_ratio_lst, Size=ais_size_status_lst)    
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 
    ais_size_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Low"], long_ratio_df$Ratio[long_ratio_df$value=="High"])
    ais_size_pvals <- append(ais_size_pvals, ais_size_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="High"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Low"]))
        ais_size_cmps <- append(ais_size_cmps, TRUE)
    else
        ais_size_cmps <- append(ais_size_cmps, FALSE)  
    ## MIA Size
    cell_ratio_df <- data.frame(Ratio=mia_size_ratio_lst, Size=mia_size_status_lst)    
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 
    mia_size_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Low"], long_ratio_df$Ratio[long_ratio_df$value=="High"])
    mia_size_pvals <- append(mia_size_pvals, mia_size_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="High"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Low"]))
        mia_size_cmps <- append(mia_size_cmps, TRUE)
    else
        mia_size_cmps <- append(mia_size_cmps, FALSE) 
    ## ADC Size
    cell_ratio_df <- data.frame(Ratio=adc_size_ratio_lst, Size=adc_size_status_lst)    
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 
    adc_size_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Low"], long_ratio_df$Ratio[long_ratio_df$value=="High"])
    adc_size_pvals <- append(adc_size_pvals, adc_size_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="High"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Low"]))
        adc_size_cmps <- append(adc_size_cmps, TRUE)
    else
        adc_size_cmps <- append(adc_size_cmps, FALSE)      
}

p_val_df <- data.frame(Age=age_pvals, Gender=gender_pvals, Race=race_pvals, Smoke=smoke_pvals,
                       Normal_AAH=normal_aah_pvals, Normal_AIS=normal_ais_pvals, Normal_MIA=normal_mia_pvals, Normal_ADC=normal_adc_pvals,
                       AAH_AIS=aah_ais_pvals, AAH_MIA=aah_mia_pvals, AAH_ADC=aah_adc_pvals, AIS_MIA=ais_mia_pvals,
                       AIS_ADC=ais_adc_pvals, MIA_ADC=mia_adc_pvals, AAH_TMB=aah_tmb_pvals, AIS_TMB=ais_tmb_pvals, 
                       MIA_TMB=mia_tmb_pvals, ADC_TMB=adc_tmb_pvals, MIA_EGFR_KRAS=mia_egfr_kras_pvals, 
                       MIA_EGFR_OTHER=mia_egfr_other_pvals, MIA_KRAS_OTHER=mia_kras_other_pvals, 
                       ADC_EGFR_KRAS=adc_egfr_kras_pvals, ADC_EGFR_OTHER=adc_egfr_other_pvals, ADC_KRAS_OTHER=adc_kras_other_pvals, 
                       AAH_Size=aah_size_pvals, AIS_Size=ais_size_pvals, MIA_Size=mia_size_pvals, ADC_Size=adc_size_pvals, 
                       row.names = all_cell_lst)
p_cmp_df <- data.frame(Age=age_cmps, Gender=gender_cmps, Race=race_cmps, Smoke=smoke_cmps,
                       Normal_AAH=normal_aah_cmps, Normal_AIS=normal_ais_cmps, Normal_MIA=normal_mia_cmps, Normal_ADC=normal_adc_cmps,
                       AAH_AIS=aah_ais_cmps, AAH_MIA=aah_mia_cmps, AAH_ADC=aah_adc_cmps, AIS_MIA=ais_mia_cmps,
                       AIS_ADC=ais_adc_cmps, MIA_ADC=mia_adc_cmps, AAH_TMB=aah_tmb_cmps, AIS_TMB=ais_tmb_cmps,
                       MIA_TMB=mia_tmb_cmps, ADC_TMB=adc_tmb_cmps, MIA_EGFR_KRAS=mia_egfr_kras_cmps, 
                       MIA_EGFR_OTHER=mia_egfr_other_cmps, MIA_KRAS_OTHER=mia_kras_other_cmps, 
                       ADC_EGFR_KRAS=adc_egfr_kras_cmps, ADC_EGFR_OTHER=adc_egfr_other_cmps, ADC_KRAS_OTHER=adc_kras_other_cmps, 
                       AAH_Size=aah_size_cmps, AIS_Size=ais_size_cmps, MIA_Size=mia_size_cmps, ADC_Size=adc_size_cmps, row.names = all_cell_lst)
var_order <- c("Age", "Gender", "Race", "Smoke", "Normal_AAH", "Normal_AIS", "Normal_MIA", "Normal_ADC", "AAH_AIS", "AAH_MIA", 
               "AAH_ADC", "AIS_MIA", "AIS_ADC", "MIA_ADC", "AAH_TMB", "AIS_TMB", "MIA_TMB", "ADC_TMB", "MIA_EGFR_KRAS", 
               "MIA_EGFR_OTHER", "MIA_KRAS_OTHER", "ADC_EGFR_KRAS", "ADC_EGFR_OTHER", "ADC_KRAS_OTHER",
               "AAH_Size", "AIS_Size", "MIA_Size", "ADC_Size")

group_df <- p_val_df %>% rownames_to_column(var = "CellType") %>% gather(Vars, Pvals, -CellType) 
group_df$AdjustPs <- p.adjust(group_df$Pvals, method = "fdr")
adjust_df <- group_df %>% mutate(gGroup = case_when(AdjustPs > 0.05 ~ 'NS', AdjustPs > 0.01 ~ '*', .default = "**"))
adjust_df$gGroup <- factor(adjust_df$gGroup, levels=c("NS", "*", "**"))
cmp_df <- p_cmp_df %>% rownames_to_column(var = "CellType") %>% gather(Vars, Cmps, -CellType) 
adjust_df$Cmps <- as.factor(cmp_df$Cmps) 

adjust_df %>% ggplot() +
    geom_point(aes(x = factor(Vars, level=var_order), y = factor(CellType, level=rev(all_cell_lst)),
                   size = p.to.Z(Pvals), shape = gGroup, color = Cmps)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))



