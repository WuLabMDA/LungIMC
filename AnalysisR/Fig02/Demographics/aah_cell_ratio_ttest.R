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
filter_df <- roi_df[roi_df$ROI_Diag %in% c("Normal", "AAH"),]
filter_df$Core_Margin <- filter_df$ROI_Location
filter_df$Core_Margin <- replace(filter_df$Core_Margin, filter_df$Core_Margin == "Tumor", "AAH")
filter_df$Core_Margin <- replace(filter_df$Core_Margin, filter_df$Core_Margin == "Normal", "DistantNormal")

all_cell_lst <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell",
                  "CD8-T-Cell", "CD4-T-Cell",  "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte",
                  "MDSC", "Fibroblast", "Undefined")
# interested ROI lst
roi_lst <- filter_df$ROI_ID
cell_id_lst <- rownames(colData(spe))
celltype_lst <- spe$celltype


dist_adja_pvals <- c()
dist_adja_cmps <- c()
dist_aah_pvals <- c()
dist_aah_cmps <- c()
adja_aah_pvals <- c()
adja_aah_cmps <- c()


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
    
    # Distant vs AAH
    dist_aah_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"], long_ratio_df$Ratio[long_ratio_df$value=="AAH"])
    dist_aah_pvals <- append(dist_aah_pvals, dist_aah_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="AAH"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="DistantNormal"]))
        dist_aah_cmps <- append(dist_aah_cmps, TRUE)
    else
        dist_aah_cmps <- append(dist_aah_cmps, FALSE)      
    
    # Adjacent vs AAH
    adja_aah_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"], long_ratio_df$Ratio[long_ratio_df$value=="AAH"])
    adja_aah_pvals <- append(adja_aah_pvals, adja_aah_ttest$p.value)
    if (mean(long_ratio_df$Ratio[long_ratio_df$value=="AAH"]) >= mean(long_ratio_df$Ratio[long_ratio_df$value=="AdjacentNormal"]))
        adja_aah_cmps <- append(adja_aah_cmps, TRUE)
    else
        adja_aah_cmps <- append(adja_aah_cmps, FALSE)         
}


p_val_df <- data.frame(Adja_Dist=dist_adja_pvals, AAH_Adja=adja_aah_pvals, row.names = all_cell_lst)
p_cmp_df <- data.frame(Adja_Dist=dist_adja_cmps, AAH_Adja=adja_aah_cmps, row.names = all_cell_lst)
var_order <- c("Adja_Dist", "AAH_Adja")

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