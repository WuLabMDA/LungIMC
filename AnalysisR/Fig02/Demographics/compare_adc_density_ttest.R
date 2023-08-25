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

# load cell phenotype information
spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

# load ROI information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_name <- "ADC_ROI_Info.xlsx"
roi_info_path <- file.path(metadata_dir, roi_info_name)
adc_df <- read_excel(roi_info_path)


all_cell_lst <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell",
                  "CD8-T-Cell", "CD4-T-Cell",  "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte",
                  "MDSC", "Fibroblast", "Undefined")
cell_id_lst <- rownames(colData(spe))
celltype_lst <- spe$celltype

# ADC
adc_dist_adja_pvals <- c()
adc_dist_adja_cmps <- c()
adc_dist_marg_pvals <- c()
adc_dist_marg_cmps <- c()
adc_dist_core_pvals <- c()
adc_dist_core_cmps <- c()
adc_adja_marg_pvals <- c()
adc_adja_marg_cmps <- c()
adc_adja_core_pvals <- c()
adc_adja_core_cmps <- c()
adc_marg_core_pvals <- c()
adc_marg_core_cmps <- c()

# ADC refining
normal_indices <- adc_df$ROI_Location %in% c("AdjacentNormal", "DistantNormal")
adc_df$Core_Margin[normal_indices] <- adc_df$ROI_Location[normal_indices]
adc_roi_lst <- adc_df$ROI_ID
roi_areas <- adc_df$Area

for (cell_type in all_cell_lst) {
    # ADC
    ratio_lst <- c()
    location_lst <- c()
    for (ir in 1:length(adc_roi_lst)) {
        cur_roi_cell_indices <- which(startsWith(cell_id_lst, adc_roi_lst[ir]))
        cur_roi_cell_types <- celltype_lst[cur_roi_cell_indices]
        roi_cell_ratio <- sum(cur_roi_cell_types == cell_type) * 1000000 / roi_areas[ir]
        ratio_lst <- append(ratio_lst, roi_cell_ratio)
        location_lst <- append(location_lst, adc_df$Core_Margin[ir])
    }
    
    
    cell_ratio_df <- data.frame(Ratio=ratio_lst, Loc=location_lst)
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio)
    # Distant vs Adjacent
    dist_adja_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"])
    adc_dist_adja_pvals <- append(adc_dist_adja_pvals, dist_adja_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        adc_dist_adja_cmps <- append(adc_dist_adja_cmps, TRUE)
    else
        adc_dist_adja_cmps <- append(adc_dist_adja_cmps, FALSE)
    # Distant vs Margin
    dist_marg_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Margin"])
    adc_dist_marg_pvals <- append(adc_dist_marg_pvals, dist_marg_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Margin"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        adc_dist_marg_cmps <- append(adc_dist_marg_cmps, TRUE)
    else
        adc_dist_marg_cmps <- append(adc_dist_marg_cmps, FALSE)
    # Distant vs Core
    dist_core_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Core"])
    adc_dist_core_pvals <- append(adc_dist_core_pvals, dist_core_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Core"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        adc_dist_core_cmps <- append(adc_dist_core_cmps, TRUE)
    else
        adc_dist_core_cmps <- append(adc_dist_core_cmps, FALSE)
    # Adjacent vs Margin
    adja_marg_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Margin"])
    adc_adja_marg_pvals <- append(adc_adja_marg_pvals, adja_marg_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Margin"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]))
        adc_adja_marg_cmps <- append(adc_adja_marg_cmps, TRUE)
    else
        adc_adja_marg_cmps <- append(adc_adja_marg_cmps, FALSE)
    # Adjacent vs Core
    adja_core_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Core"])
    adc_adja_core_pvals <- append(adc_adja_core_pvals, adja_core_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Core"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]))
        adc_adja_core_cmps <- append(adc_adja_core_cmps, TRUE)
    else
        adc_adja_core_cmps <- append(adc_adja_core_cmps, FALSE)
    # Margin vs Core
    marg_core_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Margin"], long_ratio_df$Ratio[long_ratio_df$value=="Core"])
    adc_marg_core_pvals <- append(adc_marg_core_pvals, marg_core_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Core"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Margin"]))
        adc_marg_core_cmps <- append(adc_marg_core_cmps, TRUE)
    else
        adc_marg_core_cmps <- append(adc_marg_core_cmps, FALSE)
}


p_val_df <- data.frame(ADC_Adja_Dist=adc_dist_adja_pvals, ADC_Marg_Adja=adc_adja_marg_pvals, ADC_Core_Marg=adc_marg_core_pvals, row.names = all_cell_lst)
p_cmp_df <- data.frame(ADC_Adja_Dist=adc_dist_adja_cmps, ADC_Marg_Adja=adc_adja_marg_cmps, ADC_Core_Marg=adc_marg_core_cmps, row.names = all_cell_lst)
var_order <- c("ADC_Adja_Dist", "ADC_Marg_Adja", "ADC_Core_Marg")

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