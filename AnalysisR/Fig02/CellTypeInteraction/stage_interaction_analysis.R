library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
cell_spatial_path <- file.path(phenotype_dir, paste0("InteractionsTestIter200", ".RData"))
load(cell_spatial_path)

## load ROI diagnosis information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_path <- file.path(metadata_dir, "ROI_Info.xlsx")
roi_meta_info <- read.xlsx(roi_info_path)

# Normal/AAH/AIS/MIA/ADC
path_stage <- "Normal"
if (path_stage == "Normal") {
    subset_roi_info <- subset(roi_meta_info, ROI_Diag==path_stage)
} else {
    subset_roi_info <- subset(roi_meta_info, ROI_Diag==path_stage & ROI_Location=="Tumor")
}
subset_roi_lst <- subset_roi_info$ROI_ID
subset_out <- interaction_out[interaction_out$group_by %in% subset_roi_lst, ]

from_order <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell", "CD8-T-Cell", "CD4-T-Cell", 
                "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte", "MDSC", "Fibroblast", "Undefined")
to_order <- c("Undefined", "Fibroblast", "MDSC", "Monocyte", "Macrophage", "Proliferating-Cell", "T-Reg-Cell", "CD4-T-Cell", "CD8-T-Cell", 
              "Endothelial-Cell", "Dendritic-Cell", "NK-Cell", "Neutrophil", "B-Cell", "Epithelial-Cell")

max_per_val <- 1.000
min_per_val <- -0.428

subset_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
    summarize(per_sigval = sum(sigval, na.rm = TRUE)/length(subset_roi_lst), sum_sigval = sum(sigval, na.rm = TRUE)) %>%
    mutate(across(starts_with("per"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val))) %>%
    mutate(from_label=factor(from_label, levels=from_order)) %>%
    mutate(to_label=factor(to_label, levels=to_order)) %>%
    ggplot() +
    geom_point(aes(from_label, to_label, colour=per_sigval, size=sum_sigval)) +
    scale_color_gradient2(low = "blue2", mid = "white", high = "red2") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = path_stage)
