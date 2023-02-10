library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
cell_spatial_path <- file.path(phenotype_dir, paste0("InteractionAnalysisIter100", ".RData"))
load(cell_spatial_path)

## load ROI diangoiss information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_path <- file.path(metadata_dir, "ROI_Info.xlsx")
roi_meta_info <- read.xlsx(roi_info_path)

# # filter AAH/AIS/MIA/ADC
# subset_roi_info <- subset(roi_meta_info, ROI_Diag=="AAH" & ROI_Location=="Tumor")
# Normal
subset_roi_info <- subset(roi_meta_info, ROI_Diag=="Normal")

subset_roi_lst <- subset_roi_info$ROI_ID
subset_out <- out[out$group_by %in% subset_roi_lst, ]
subset_roi_num <- length(subset_roi_lst)

celltype_order <- c("Monocytes", "Unknown", "Macrophage", "Endothelial cells", "Epithelial",
                    "Neutrophils", "Smooth muscle/Stromal", "CD4 T cell", "CD8 T cells", "Dendritic cell",
                    "MDSC", "NK cell", "T-reg cells", "B cells")

subset_out %>% as_tibble() %>%
    group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / subset_roi_num) %>%
    mutate(from_label=factor(from_label, levels=celltype_order)) %>% 
    mutate(to_label=factor(to_label, levels=celltype_order)) %>% 
    ggplot() +
    geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
    scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

