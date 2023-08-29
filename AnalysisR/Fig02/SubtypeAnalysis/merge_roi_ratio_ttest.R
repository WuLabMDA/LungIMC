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
spe_celltype_name <-"lung_spe_32_cell_subtypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)


# load ROI information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_name <- "ROI_Info_Aggregation.csv"
roi_info_path <- file.path(metadata_dir, roi_info_name)
roi_df <- read_csv(roi_info_path)


# list all cell types
all_cell_lst <- c("Ki67+ Epithelial", "Ki67+ PDL1+ Epithelial", "CD73+ Epithelial", "Other Epithelial", 
                  "Ki67+ B-Cells", "Ki67- B-Cells", "Ki67+ NK Cells", "Ki67- NK Cells",
                  "Ki67+ Dendritic Cells", "HLADR+ Dendritic Cells", "Other Dendritic Cells", "Endothelial-Cell",
                  "Cytotoxic CD8 T-Cells", "Memory CD8 T-Cells", "Exhausted CD8 T-Cells", "Ki67+ CD8 T-Cells", "Naive CD8 T-Cells",
                  "Memory CD4 T-Cells", "Exhausted CD4 T-Cells", "Ki67+ CD4 T-Cells", "Naive CD4 T-Cells",
                  "Ki67+ Treg-Cells", "Ki67- Treg-Cells", "Proliferating-Cell", 
                  "CD163+ Macrophages", "Ki67+ Macrophages", "CD163- Macrophages",
                  "Neutrophil", "Monocyte", "MDSC", "Fibroblast", "Undefined")
cell_id_lst <- rownames(colData(spe))
celltype_lst <- spe$cellsubtype

# AAH
aah_dist_adja_pvals <- c()
aah_dist_adja_cmps <- c()
aah_dist_aah_pvals <- c()
aah_dist_aah_cmps <- c()
aah_adja_aah_pvals <- c()
aah_adja_aah_cmps <- c()
# AIS
ais_dist_adja_pvals <- c()
ais_dist_adja_cmps <- c()
ais_dist_marg_pvals <- c()
ais_dist_marg_cmps <- c()
ais_dist_core_pvals <- c()
ais_dist_core_cmps <- c()
ais_adja_marg_pvals <- c()
ais_adja_marg_cmps <- c()
ais_adja_core_pvals <- c()
ais_adja_core_cmps <- c()
ais_marg_core_pvals <- c()
ais_marg_core_cmps <- c()
# MIA
mia_dist_adja_pvals <- c()
mia_dist_adja_cmps <- c()
mia_dist_marg_pvals <- c()
mia_dist_marg_cmps <- c()
mia_dist_core_pvals <- c()
mia_dist_core_cmps <- c()
mia_adja_marg_pvals <- c()
mia_adja_marg_cmps <- c()
mia_adja_core_pvals <- c()
mia_adja_core_cmps <- c()
mia_marg_core_pvals <- c()
mia_marg_core_cmps <- c()
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


# AAH filter
aah_filter_df <- roi_df[roi_df$ROI_Diag %in% c("Normal", "AAH"),]
aah_filter_df$Core_Margin <- aah_filter_df$ROI_Location
aah_filter_df$Core_Margin <- replace(aah_filter_df$Core_Margin, aah_filter_df$Core_Margin == "Tumor", "AAH")
aah_filter_df$Core_Margin <- replace(aah_filter_df$Core_Margin, aah_filter_df$Core_Margin == "Normal", "DistantNormal")
# AIS filter
ais_df <- roi_df[roi_df$ROI_Diag=="AIS",]
normal_indices <- ais_df$ROI_Location %in% c("AdjacentNormal", "DistantNormal")
ais_df$Core_Margin[normal_indices] <- ais_df$ROI_Location[normal_indices]
ais_filter_df <- ais_df[!is.na(ais_df$Core_Margin),]
# MIA filter
mia_df <- roi_df[roi_df$ROI_Diag=="MIA",]
normal_indices <- mia_df$ROI_Location %in% c("AdjacentNormal", "DistantNormal")
mia_df$Core_Margin[normal_indices] <- mia_df$ROI_Location[normal_indices]
mia_filter_df <- mia_df[!is.na(mia_df$Core_Margin),]
# ADC filter
adc_df <- roi_df[roi_df$ROI_Diag=="ADC",]
normal_indices <- adc_df$ROI_Location %in% c("AdjacentNormal", "DistantNormal")
adc_df$Core_Margin[normal_indices] <- adc_df$ROI_Location[normal_indices]
adc_filter_df <- adc_df[!is.na(adc_df$Core_Margin),]


# interested ROI lst
aah_roi_lst <- aah_filter_df$ROI_ID
ais_roi_lst <- ais_filter_df$ROI_ID
mia_roi_lst <- mia_filter_df$ROI_ID
adc_roi_lst <- adc_filter_df$ROI_ID


for (cell_type in all_cell_lst) {
    # AAH
    ratio_lst <- c()
    location_lst <- c()
    for (ir in 1:length(aah_roi_lst)) {
        cur_roi_cell_indices <- which(startsWith(cell_id_lst, aah_roi_lst[ir]))
        cur_roi_cell_types <- celltype_lst[cur_roi_cell_indices]
        roi_cell_ratio <- sum(cur_roi_cell_types == cell_type) / length(cur_roi_cell_types)
        ratio_lst <- append(ratio_lst, roi_cell_ratio)
        location_lst <- append(location_lst, aah_filter_df$Core_Margin[ir])
    }
    
    cell_ratio_df <- data.frame(Ratio=ratio_lst, Loc=location_lst)
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 
    
    # Distant vs Adjacent
    dist_adja_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"])
    aah_dist_adja_pvals <- append(aah_dist_adja_pvals, dist_adja_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        aah_dist_adja_cmps <- append(aah_dist_adja_cmps, TRUE)
    else
        aah_dist_adja_cmps <- append(aah_dist_adja_cmps, FALSE)   
    # Distant vs AAH
    dist_aah_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="AAH"])
    aah_dist_aah_pvals <- append(aah_dist_aah_pvals, dist_aah_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="AAH"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        aah_dist_aah_cmps <- append(aah_dist_aah_cmps, TRUE)
    else
        aah_dist_aah_cmps <- append(aah_dist_aah_cmps, FALSE)      
    # Adjacent vs AAH
    adja_aah_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"], long_ratio_df$Ratio[long_ratio_df$value=="AAH"])
    aah_adja_aah_pvals <- append(aah_adja_aah_pvals, adja_aah_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="AAH"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]))
        aah_adja_aah_cmps <- append(aah_adja_aah_cmps, TRUE)
    else
        aah_adja_aah_cmps <- append(aah_adja_aah_cmps, FALSE)    
    
    # AIS
    ratio_lst <- c()
    location_lst <- c()
    for (ir in 1:length(ais_roi_lst)) {
        cur_roi_cell_indices <- which(startsWith(cell_id_lst, ais_roi_lst[ir]))
        cur_roi_cell_types <- celltype_lst[cur_roi_cell_indices]
        roi_cell_ratio <- sum(cur_roi_cell_types == cell_type) / length(cur_roi_cell_types)
        ratio_lst <- append(ratio_lst, roi_cell_ratio)
        location_lst <- append(location_lst, ais_filter_df$Core_Margin[ir])
    }
    
    cell_ratio_df <- data.frame(Ratio=ratio_lst, Loc=location_lst)
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio)
    # Distant vs Adjacent
    dist_adja_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"])
    ais_dist_adja_pvals <- append(ais_dist_adja_pvals, dist_adja_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        ais_dist_adja_cmps <- append(ais_dist_adja_cmps, TRUE)
    else
        ais_dist_adja_cmps <- append(ais_dist_adja_cmps, FALSE)
    # Distant vs Margin
    dist_marg_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Margin"])
    ais_dist_marg_pvals <- append(ais_dist_marg_pvals, dist_marg_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Margin"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        ais_dist_marg_cmps <- append(ais_dist_marg_cmps, TRUE)
    else
        ais_dist_marg_cmps <- append(ais_dist_marg_cmps, FALSE)
    # Distant vs Core
    dist_core_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Core"])
    ais_dist_core_pvals <- append(ais_dist_core_pvals, dist_core_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Core"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        ais_dist_core_cmps <- append(ais_dist_core_cmps, TRUE)
    else
        ais_dist_core_cmps <- append(ais_dist_core_cmps, FALSE)
    # Adjacent vs Margin
    adja_marg_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Margin"])
    ais_adja_marg_pvals <- append(ais_adja_marg_pvals, adja_marg_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Margin"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]))
        ais_adja_marg_cmps <- append(ais_adja_marg_cmps, TRUE)
    else
        ais_adja_marg_cmps <- append(ais_adja_marg_cmps, FALSE)
    # Adjacent vs Core
    adja_core_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Core"])
    ais_adja_core_pvals <- append(ais_adja_core_pvals, adja_core_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Core"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]))
        ais_adja_core_cmps <- append(ais_adja_core_cmps, TRUE)
    else
        ais_adja_core_cmps <- append(ais_adja_core_cmps, FALSE)
    # Margin vs Core
    marg_core_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Margin"], long_ratio_df$Ratio[long_ratio_df$value=="Core"])
    ais_marg_core_pvals <- append(ais_marg_core_pvals, marg_core_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Core"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Margin"]))
        ais_marg_core_cmps <- append(ais_marg_core_cmps, TRUE)
    else
        ais_marg_core_cmps <- append(ais_marg_core_cmps, FALSE)
    
    # MIA
    ratio_lst <- c()
    location_lst <- c()
    for (ir in 1:length(mia_roi_lst)) {
        cur_roi_cell_indices <- which(startsWith(cell_id_lst, mia_roi_lst[ir]))
        cur_roi_cell_types <- celltype_lst[cur_roi_cell_indices]
        roi_cell_ratio <- sum(cur_roi_cell_types == cell_type) / length(cur_roi_cell_types)
        ratio_lst <- append(ratio_lst, roi_cell_ratio)
        location_lst <- append(location_lst, mia_filter_df$Core_Margin[ir])
    }
    
    cell_ratio_df <- data.frame(Ratio=ratio_lst, Loc=location_lst)
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio)
    # Distant vs Adjacent
    dist_adja_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"])
    mia_dist_adja_pvals <- append(mia_dist_adja_pvals, dist_adja_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        mia_dist_adja_cmps <- append(mia_dist_adja_cmps, TRUE)
    else
        mia_dist_adja_cmps <- append(mia_dist_adja_cmps, FALSE)
    # Distant vs Margin
    dist_marg_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Margin"])
    mia_dist_marg_pvals <- append(mia_dist_marg_pvals, dist_marg_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Margin"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        mia_dist_marg_cmps <- append(mia_dist_marg_cmps, TRUE)
    else
        mia_dist_marg_cmps <- append(mia_dist_marg_cmps, FALSE)
    # Distant vs Core
    dist_core_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Core"])
    mia_dist_core_pvals <- append(mia_dist_core_pvals, dist_core_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Core"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        mia_dist_core_cmps <- append(mia_dist_core_cmps, TRUE)
    else
        mia_dist_core_cmps <- append(mia_dist_core_cmps, FALSE)
    # Adjacent vs Margin
    adja_marg_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Margin"])
    mia_adja_marg_pvals <- append(mia_adja_marg_pvals, adja_marg_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Margin"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]))
        mia_adja_marg_cmps <- append(mia_adja_marg_cmps, TRUE)
    else
        mia_adja_marg_cmps <- append(mia_adja_marg_cmps, FALSE)
    # Adjacent vs Core
    adja_core_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Core"])
    mia_adja_core_pvals <- append(mia_adja_core_pvals, adja_core_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Core"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]))
        mia_adja_core_cmps <- append(mia_adja_core_cmps, TRUE)
    else
        mia_adja_core_cmps <- append(mia_adja_core_cmps, FALSE)
    # Margin vs Core
    marg_core_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Margin"], long_ratio_df$Ratio[long_ratio_df$value=="Core"])
    mia_marg_core_pvals <- append(mia_marg_core_pvals, marg_core_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Core"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Margin"]))
        mia_marg_core_cmps <- append(mia_marg_core_cmps, TRUE)
    else
        mia_marg_core_cmps <- append(mia_marg_core_cmps, FALSE)
    # ADC
    ratio_lst <- c()
    location_lst <- c()
    for (ir in 1:length(adc_roi_lst)) {
        cur_roi_cell_indices <- which(startsWith(cell_id_lst, adc_roi_lst[ir]))
        cur_roi_cell_types <- celltype_lst[cur_roi_cell_indices]
        roi_cell_ratio <- sum(cur_roi_cell_types == cell_type) / length(cur_roi_cell_types)
        ratio_lst <- append(ratio_lst, roi_cell_ratio)
        location_lst <- append(location_lst, adc_filter_df$Core_Margin[ir])
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


p_val_df <- data.frame(AAH_Adja_Dist=aah_dist_adja_pvals, AAH_AAH_Adja=aah_adja_aah_pvals, 
                       AIS_Adja_Dist=ais_dist_adja_pvals, AIS_Marg_Adja=ais_adja_marg_pvals, AIS_Core_Marg=ais_marg_core_pvals, 
                       MIA_Adja_Dist=mia_dist_adja_pvals, MIA_Marg_Adja=mia_adja_marg_pvals, MIA_Core_Marg=mia_marg_core_pvals,
                       ADC_Adja_Dist=adc_dist_adja_pvals, ADC_Marg_Adja=adc_adja_marg_pvals, ADC_Core_Marg=adc_marg_core_pvals,
                       row.names = all_cell_lst)
p_cmp_df <- data.frame(AAH_Adja_Dist=aah_dist_adja_cmps, AAH_AAH_Adja=aah_adja_aah_cmps, 
                       AIS_Adja_Dist=ais_dist_adja_cmps, AIS_Marg_Adja=ais_adja_marg_cmps, AIS_Core_Marg=ais_marg_core_cmps,
                       MIA_Adja_Dist=mia_dist_adja_cmps, MIA_Marg_Adja=mia_adja_marg_cmps, MIA_Core_Marg=mia_marg_core_cmps,
                       ADC_Adja_Dist=adc_dist_adja_cmps, ADC_Marg_Adja=adc_adja_marg_cmps, ADC_Core_Marg=adc_marg_core_cmps,
                       row.names = all_cell_lst)
var_order <- c("AAH_Adja_Dist", "AAH_AAH_Adja", "AIS_Adja_Dist", "AIS_Marg_Adja", "AIS_Core_Marg",
               "MIA_Adja_Dist", "MIA_Marg_Adja", "MIA_Core_Marg", "ADC_Adja_Dist", "ADC_Marg_Adja", "ADC_Core_Marg")


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
