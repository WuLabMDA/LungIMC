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
# load slide mutation information
slide_mutation_path <- file.path(metadata_dir, "TMB", "LungSlideMutation.csv")
slide_mutation_df <- read_csv(slide_mutation_path)
# load ROI information
roi_info_name <- "ROI_Info"
roi_info_path <- file.path(metadata_dir, paste0(roi_info_name, ".xlsx"))
roi_df <- read_excel(roi_info_path)
# load lesion information
slide_info_path <- file.path(metadata_dir, "Lesion_Info.xlsx")
slide_info_df <- read_excel(slide_info_path)


# pathological stages
interested_stage <- "MIA"

## filtering rois
roi_lst <- roi_df$ROI_ID
stage_roi_vec <- c()
for (ind in 1:nrow(roi_df)) {
    roi_location <- roi_df$ROI_Location[ind]
    roi_diag <- roi_df$ROI_Diag[ind]
    if (roi_location == "Tumor" & roi_diag == interested_stage)
        stage_roi_vec <- append(stage_roi_vec, roi_df$ROI_ID[ind])
}

## obtain cells inside tumor lesions
lesion_indices = c()
cell_id_lst <- rownames(colData(spe))
for (cur_roi in stage_roi_vec) 
    lesion_indices <- append(lesion_indices, which(startsWith(cell_id_lst, cur_roi)))
cell_id_lst <- cell_id_lst[lesion_indices]
celltype_lst <- spe$celltype
celltype_lst <- celltype_lst[lesion_indices]

# filter slide
mutation_diag <- c()
for (slide_name in slide_mutation_df$Slides) {
    filter_slide <- slide_info_df[slide_info_df$Slide_ID == slide_name, ]
    if (nrow(filter_slide) != 1) {
        print(slide_name)
        print(nrow(filter_slide))
    }
    mutation_diag <- append(mutation_diag, filter_slide$Slide_Diag)
}
slide_mutation_df$Diag <- mutation_diag
stage_mutation_df <- slide_mutation_df[slide_mutation_df$Diag==interested_stage, ]


all_cell_lst <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell",
                  "CD8-T-Cell", "CD4-T-Cell",  "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte",
                  "MDSC", "Fibroblast", "Undefined")

egfr_kras_pvals <- c()
egfr_kras_cmps <- c()
egfr_other_pvals <- c()
egfr_other_cmps <- c()
kras_other_pvals <- c()
kras_other_cmps <- c()
# tp53_other_pvals <- c()
# tp53_other_cmps <- c()

for (cell_type in all_cell_lst) {
    # collect information
    ratio_lst <- c()
    mutation_lst <- c()    
    # iterate by slide
    stage_slide_lst <- stage_mutation_df$Slides
    for (ir in 1:length(stage_slide_lst)) {
        cur_slide <- stage_slide_lst[ir]
        slide_cell_indices <- which(startsWith(cell_id_lst, cur_slide))
        slide_celltypes <- celltype_lst[slide_cell_indices]
        slide_cell_ratio <- sum(slide_celltypes == cell_type) / length(slide_celltypes)
        ratio_lst <- append(ratio_lst, slide_cell_ratio)
    }    

    cell_ratio_df <- data.frame(Ratio=ratio_lst, EGFR_KRAS=stage_mutation_df$`EGFR-KRAS`)    
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 

    egfr_kras_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="EGFR"], long_ratio_df$Ratio[long_ratio_df$value=="KRAS"])
    egfr_kras_pvals <- append(egfr_kras_pvals, egfr_kras_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="EGFR"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="KRAS"]))
        egfr_kras_cmps <- append(egfr_kras_cmps, TRUE)
    else
        egfr_kras_cmps <- append(egfr_kras_cmps, FALSE)      
    
    egfr_other_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="EGFR"], long_ratio_df$Ratio[long_ratio_df$value=="Other"])
    egfr_other_pvals <- append(egfr_other_pvals, egfr_other_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="EGFR"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Other"]))
        egfr_other_cmps <- append(egfr_other_cmps, TRUE)
    else
        egfr_other_cmps <- append(egfr_other_cmps, FALSE)  
    
    kras_other_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="KRAS"], long_ratio_df$Ratio[long_ratio_df$value=="Other"])
    kras_other_pvals <- append(kras_other_pvals, kras_other_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="KRAS"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Other"]))
        kras_other_cmps <- append(kras_other_cmps, TRUE)
    else
        kras_other_cmps <- append(kras_other_cmps, FALSE)           
}

p_val_df <- data.frame(EGFR_KRAS=egfr_kras_pvals, EGFR_Other=egfr_other_pvals, KRAS_Other=kras_other_pvals, row.names = all_cell_lst)
p_cmp_df <- data.frame(EGFR_KRAS=egfr_kras_cmps, EGFR_Other=egfr_other_cmps, KRAS_Other=kras_other_cmps, row.names = all_cell_lst)
var_order <- c("EGFR_KRAS", "EGFR_Other", "KRAS_Other")

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

# adjust_df %>% ggplot() +
#     geom_point(aes(x = factor(Vars, level=var_order), y = factor(CellType, level=rev(all_cell_lst)),
#                    size = p.to.Z(Pvals), color = Cmps)) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))

