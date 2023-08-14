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
        cur_size <- slide_info_df$TumorSize[ir]
        cur_diag <- slide_info_df$Slide_Diag[ir]
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

p_val_df <- data.frame(AAH_Size=aah_size_pvals, AIS_Size=ais_size_pvals, 
                       MIA_Size=mia_size_pvals, ADC_Size=adc_size_pvals, row.names = all_cell_lst)
p_cmp_df <- data.frame(AAH_Size=aah_size_cmps, AIS_Size=ais_size_cmps,
                       MIA_Size=mia_size_cmps, ADC_Size=adc_size_cmps, row.names = all_cell_lst)
var_order <- c("AAH_Size", "AIS_Size", "MIA_Size", "ADC_Size")

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
