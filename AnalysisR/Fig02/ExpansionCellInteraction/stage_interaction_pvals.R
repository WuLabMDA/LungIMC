library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
celltype_expansion_dir <- file.path(phenotype_dir, "ExpansionInteraction")
threshold_val <- 50
cell_type_interaction_name <- paste0("ExpansionInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_expansion_dir, cell_type_interaction_name)
cell_spatial_path <- file.path(cell_type_interaction_path)
load(cell_spatial_path)

## load ROI diagnosis information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_path <- file.path(metadata_dir, "ROI_Info.xlsx")
roi_meta_info <- read.xlsx(roi_info_path)


from_order <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", 
                "Endothelial-Cell", "CD8-T-Cell", "CD4-T-Cell", "T-Reg-Cell", "Proliferating-Cell", 
                "Macrophage", "Monocyte", "MDSC", "Fibroblast", "Undefined")
to_order <- c("Undefined", "Fibroblast", "MDSC", "Monocyte", "Macrophage", 
              "Proliferating-Cell", "T-Reg-Cell", "CD4-T-Cell", "CD8-T-Cell", "Endothelial-Cell", 
              "Dendritic-Cell", "NK-Cell", "Neutrophil", "B-Cell", "Epithelial-Cell")

# subset roi information
normal_roi_info <- subset(roi_meta_info, ROI_Diag=="Normal")
aah_roi_info <- subset(roi_meta_info, ROI_Diag=="AAH" & ROI_Location=="Tumor")
ais_roi_info <- subset(roi_meta_info, ROI_Diag=="AIS" & ROI_Location=="Tumor")
mia_roi_info <- subset(roi_meta_info, ROI_Diag=="MIA" & ROI_Location=="Tumor")
adc_roi_info <- subset(roi_meta_info, ROI_Diag=="ADC" & ROI_Location=="Tumor")

# subset interaction
normal_subset_out <- interaction_out[interaction_out$group_by %in% normal_roi_info$ROI_ID, ]
normal_subset <- normal_subset_out %>% as_tibble() %>% group_by(from_label, to_label)
aah_subset_out <- interaction_out[interaction_out$group_by %in% aah_roi_info$ROI_ID, ]
aah_subset <- aah_subset_out %>% as_tibble() %>% group_by(from_label, to_label)
ais_subset_out <- interaction_out[interaction_out$group_by %in% ais_roi_info$ROI_ID, ]
ais_subset <- ais_subset_out %>% as_tibble() %>% group_by(from_label, to_label)
mia_subset_out <- interaction_out[interaction_out$group_by %in% mia_roi_info$ROI_ID, ]
mia_subset <- mia_subset_out %>% as_tibble() %>% group_by(from_label, to_label)
adc_subset_out <- interaction_out[interaction_out$group_by %in% adc_roi_info$ROI_ID, ]
adc_subset <- adc_subset_out %>% as_tibble() %>% group_by(from_label, to_label)

from_order <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", 
                "Endothelial-Cell", "CD8-T-Cell", "CD4-T-Cell", "T-Reg-Cell", "Proliferating-Cell", 
                "Macrophage", "Monocyte", "MDSC", "Fibroblast", "Undefined")
to_order <- c("Undefined", "Fibroblast", "MDSC", "Monocyte", "Macrophage", 
              "Proliferating-Cell", "T-Reg-Cell", "CD4-T-Cell", "CD8-T-Cell", "Endothelial-Cell", 
              "Dendritic-Cell", "NK-Cell", "Neutrophil", "B-Cell", "Epithelial-Cell")

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
        normal_pair <- normal_subset[normal_subset$from_label == p_from_label & normal_subset$to_label == p_to_label,] %>% drop_na(sigval)
        aah_pair <- aah_subset[aah_subset$from_label == p_from_label & aah_subset$to_label == p_to_label,] %>% drop_na(sigval)
        ais_pair <- ais_subset[ais_subset$from_label == p_from_label & ais_subset$to_label == p_to_label,] %>% drop_na(sigval)
        mia_pair <- mia_subset[mia_subset$from_label == p_from_label & mia_subset$to_label == p_to_label,] %>% drop_na(sigval)
        adc_pair <- adc_subset[adc_subset$from_label == p_from_label & adc_subset$to_label == p_to_label,] %>% drop_na(sigval)
        
        kw_pval <- kruskal.test(list(normal_pair$sigval, aah_pair$sigval, ais_pair$sigval, mia_pair$sigval, adc_pair$sigval))
  
        pval_df[p_to_label, p_from_label] <- kw_pval$p.value
        group_pvals[group_ind, ] <- c(p_to_label, p_from_label, kw_pval$p.value)
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
