library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)
library(NCmisc)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
celltype_expansion_dir <- file.path(phenotype_dir, "ExpansionInteraction")
threshold_val <- 24
cell_type_interaction_name <- paste0("ExpansionInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_expansion_dir, cell_type_interaction_name)
cell_spatial_path <- file.path(cell_type_interaction_path)
load(cell_spatial_path)

## load ROI diagnosis information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_path <- file.path(metadata_dir, "ROI_Info.xlsx")
roi_meta_info <- read.xlsx(roi_info_path)

## load tbm information
tmb_info_path <- file.path(metadata_dir, "TMB", "LungSlideTMB2.csv")
roi_tmb_info <- read.csv(tmb_info_path)

# AAH/AIS/MIA/ADC
path_stage <- "ADC"
subset_roi_info <- subset(roi_meta_info, ROI_Diag==path_stage & ROI_Location=="Tumor")
roi_slides <- str_extract_all(subset_roi_info$ROI_ID, ".+(?=-ROI)", simplify = TRUE)
subset_roi_info <- cbind(subset_roi_info, roi_slides)

subset_tmb_info <- filter(roi_tmb_info, Stages==path_stage)
low_subset_tmb <- filter(subset_tmb_info, TMB.Cat2=="Low")
high_subset_tmb <- filter(subset_tmb_info, TMB.Cat2=="High")


low_sub_roi_info <- filter(subset_roi_info, roi_slides %in% low_subset_tmb$Slides)
high_sub_roi_info <- filter(subset_roi_info, roi_slides %in% high_subset_tmb$Slides)
low_sub_roi_lst <- low_sub_roi_info$ROI_ID
high_sub_roi_lst <- high_sub_roi_info$ROI_ID

from_order <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", 
                "Endothelial-Cell", "CD8-T-Cell", "CD4-T-Cell", "T-Reg-Cell", "Proliferating-Cell", 
                "Macrophage", "Monocyte", "MDSC", "Fibroblast", "Undefined")
to_order <- c("Undefined", "Fibroblast", "MDSC", "Monocyte", "Macrophage", 
              "Proliferating-Cell", "T-Reg-Cell", "CD4-T-Cell", "CD8-T-Cell", "Endothelial-Cell", 
              "Dendritic-Cell", "NK-Cell", "Neutrophil", "B-Cell", "Epithelial-Cell")

low_subset_out <- interaction_out[interaction_out$group_by %in% low_sub_roi_lst, ]
high_subset_out <- interaction_out[interaction_out$group_by %in% high_sub_roi_lst, ]


low_subset <- low_subset_out %>% as_tibble() %>% group_by(from_label, to_label)
high_subset <- high_subset_out %>% as_tibble() %>% group_by(from_label, to_label)


pval_df <- data.frame(matrix(nrow = length(from_order), ncol = length(from_order)))
rownames(pval_df) <- from_order
colnames(pval_df) <- from_order

group_pvals <- data.frame(matrix(nrow = length(from_order) * length(from_order), ncol = 3))
colnames(group_pvals) <- c("FromType", "ToType", "pVal")

p_from_order <- from_order
p_to_order <- from_order

group_ind <- 1
for (p_from_label in p_from_order) {
    for (p_to_label in p_to_order) {
        low_pair <- low_subset[low_subset$from_label == p_from_label & low_subset$to_label == p_to_label,] %>% drop_na(sigval)
        high_pair <- high_subset[high_subset$from_label == p_from_label & high_subset$to_label == p_to_label,] %>% drop_na(sigval)
        merge_sig_vals <- sigvals <- c(high_pair$sigval, low_pair$sigval)
        if (all(merge_sig_vals == merge_sig_vals[1]))
            cur_p_val <- 0.99
        else {
            test_pval <- t.test(high_pair$sigval, low_pair$sigval)
            cur_p_val <- test_pval$p.value
        }
        pval_df[p_to_label, p_from_label] <- cur_p_val
        group_pvals[group_ind, ] <- c(p_to_label, p_from_label, cur_p_val)
        group_ind <- group_ind + 1
    }
}


adjust_group <- group_pvals %>% 
    mutate(pVal = as.numeric(pVal)) %>%
    mutate(gGroup = case_when(pVal > 0.05 ~ 'NS', pVal > 0.01 ~ '*', pVal > 0.001 ~ '**', .default = "***")) %>%
    mutate(gGroup = factor(gGroup, level=c("NS", "*", "**", "***")))


adjust_group %>% ggplot() +
    geom_point(aes(x = factor(ToType, level=from_order), y = factor(FromType, level=to_order), 
                   size = p.to.Z(pVal), col = gGroup)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 