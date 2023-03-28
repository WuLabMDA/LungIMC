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
roi_info_name <- "ROI_Info"
roi_info_path <- file.path(metadata_dir, paste0(roi_info_name, ".xlsx"))
roi_df <- read_excel(roi_info_path)
mia_df <- roi_df[roi_df$ROI_Diag=="MIA",]
normal_indices <- mia_df$ROI_Location %in% c("AdjacentNormal", "DistantNormal")
mia_df$Core_Margin[normal_indices] <- mia_df$ROI_Location[normal_indices]
filter_df <- mia_df[!is.na(mia_df$Core_Margin),]


all_cell_lst <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell",
                  "CD8-T-Cell", "CD4-T-Cell",  "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte",
                  "MDSC", "Fibroblast", "Undefined")
# interested ROI lst
roi_lst <- filter_df$ROI_ID
cell_id_lst <- rownames(colData(spe))
celltype_lst <- spe$celltype


dist_adja_pvals <- c()
dist_adja_cmps <- c()
dist_marg_pvals <- c()
dist_marg_cmps <- c()
dist_core_pvals <- c()
dist_core_cmps <- c()
adja_marg_pvals <- c()
adja_marg_cmps <- c()
adja_core_pvals <- c()
adja_core_cmps <- c()
marg_core_pvals <- c()
marg_core_cmps <- c()


for (cell_type in all_cell_lst) {
    ratio_lst <- c()
    location_lst <- c()
    for (ir in 1:length(roi_lst)) {
        cur_roi_cell_indices <- which(startsWith(cell_id_lst, roi_lst[ir]))
        cur_roi_cell_types <- celltype_lst[cur_roi_cell_indices]
        roi_cell_ratio <- sum(cur_roi_cell_types == cell_type) / length(cur_roi_cell_types)
        ratio_lst <- append(ratio_lst, roi_cell_ratio)
        location_lst <- append(location_lst, filter_df$Core_Margin[ir])
    }
    
    cell_ratio_df <- data.frame(Ratio=ratio_lst, Loc=location_lst)
    long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio) 
    # Distant vs Adjacent
    dist_adja_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"])
    dist_adja_pvals <- append(dist_adja_pvals, dist_adja_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        dist_adja_cmps <- append(dist_adja_cmps, TRUE)
    else
        dist_adja_cmps <- append(dist_adja_cmps, FALSE)
    # Distant vs Margin
    dist_marg_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Margin"])
    dist_marg_pvals <- append(dist_marg_pvals, dist_marg_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Margin"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        dist_marg_cmps <- append(dist_marg_cmps, TRUE)
    else
        dist_marg_cmps <- append(dist_marg_cmps, FALSE)    
    # Distant vs Core
    dist_core_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Core"])
    dist_core_pvals <- append(dist_core_pvals, dist_core_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Core"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        dist_core_cmps <- append(dist_core_cmps, TRUE)
    else
        dist_core_cmps <- append(dist_core_cmps, FALSE)
    # Adjacent vs Margin
    adja_marg_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Margin"])
    adja_marg_pvals <- append(adja_marg_pvals, adja_marg_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Margin"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]))
        adja_marg_cmps <- append(adja_marg_cmps, TRUE)
    else
        adja_marg_cmps <- append(adja_marg_cmps, FALSE)    
    # Adjacent vs Core
    adja_core_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"], long_ratio_df$Ratio[long_ratio_df$value=="Core"])
    adja_core_pvals <- append(adja_core_pvals, adja_core_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Core"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]))
        adja_core_cmps <- append(adja_core_cmps, TRUE)
    else
        adja_core_cmps <- append(adja_core_cmps, FALSE)   
    # Margin vs Core
    marg_core_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Margin"], long_ratio_df$Ratio[long_ratio_df$value=="Core"])
    marg_core_pvals <- append(marg_core_pvals, marg_core_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="Core"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="Margin"]))
        marg_core_cmps <- append(marg_core_cmps, TRUE)
    else
        marg_core_cmps <- append(marg_core_cmps, FALSE)        
}

# p_val_df <- data.frame(Adja_Dist=dist_adja_pvals, Marg_Dist=dist_marg_pvals, Core_Dist=dist_core_pvals, 
#                        Marg_Adja=adja_marg_pvals, Core_Adja=adja_core_pvals, Core_Marg=marg_core_pvals,
#                        row.names = all_cell_lst)
# p_cmp_df <- data.frame(Adja_Dist=dist_adja_cmps, Marg_Dist=dist_marg_cmps, Core_Dist=dist_core_cmps, 
#                        Marg_Adja=adja_marg_cmps, Core_Adja=adja_core_cmps, Core_Marg=marg_core_cmps,
#                        row.names = all_cell_lst)
# var_order <- c("Adja_Dist", "Marg_Dist", "Core_Dist", "Marg_Adja", "Core_Adja", "Core_Marg")

p_val_df <- data.frame(Adja_Dist=dist_adja_pvals, Marg_Adja=adja_marg_pvals, Core_Marg=marg_core_pvals, row.names = all_cell_lst)
p_cmp_df <- data.frame(Adja_Dist=dist_adja_cmps, Marg_Adja=adja_marg_cmps, Core_Marg=marg_core_cmps, row.names = all_cell_lst)
var_order <- c("Adja_Dist", "Marg_Adja", "Core_Marg")

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

