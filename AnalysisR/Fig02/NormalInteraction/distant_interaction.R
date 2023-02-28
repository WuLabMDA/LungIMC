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


# AdjacentNormal of Normal/AAH/AIS/MIA/ADC
path_stage <- "ADC"
subset_roi_info <- subset(roi_meta_info, ROI_Diag==path_stage & ROI_Location=="DistantNormal")


subset_roi_lst <- subset_roi_info$ROI_ID
subset_roi_num <- length(subset_roi_lst)
subset_out <- interaction_out[interaction_out$group_by %in% subset_roi_lst, ]

from_order <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell", "CD8-T-Cell", "CD4-T-Cell", 
                "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte", "MDSC", "Fibroblast", "Undefined")
to_order <- c("Undefined", "Fibroblast", "MDSC", "Monocyte", "Macrophage", "Proliferating-Cell", "T-Reg-Cell", "CD4-T-Cell", "CD8-T-Cell", 
              "Endothelial-Cell", "Dendritic-Cell", "NK-Cell", "Neutrophil", "B-Cell", "Epithelial-Cell")

max_per_val <- 1.000
min_per_val <- -0.428

subset_out %>% as_tibble() %>%
    group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / subset_roi_num) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val))) %>%
    mutate(from_label=factor(from_label, levels=from_order)) %>%
    mutate(to_label=factor(to_label, levels=to_order)) %>%
    ggplot() +
    geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))