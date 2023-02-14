library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
cell_spatial_path <- file.path(phenotype_dir, paste0("InteractionAnalysisIter200", ".RData"))
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

from_order <- c("Epithelial", "B cells", "Neutrophils", "NK cell", "Dendritic cell", "Endothelial cells", 
                "CD8 T cells", "CD4 T cell", "T-reg cells", "Macrophage", "Monocytes", "MDSC", 
                "Smooth muscle/Stromal", "Unknown")
to_order <- c("Unknown", "Smooth muscle/Stromal", "MDSC", "Monocytes", "Macrophage", "T-reg cells", 
              "CD4 T cell", "CD8 T cells", "Endothelial cells", "Dendritic cell", "NK cell", "Neutrophils", 
              "B cells", "Epithelial")

subset_out %>% as_tibble() %>%
    group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(subset_roi_lst)) %>%
    mutate(from_label=factor(from_label, levels=from_order)) %>% 
    mutate(to_label=factor(to_label, levels=to_order)) %>% 
    ggplot() +
    geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(title = path_stage)

