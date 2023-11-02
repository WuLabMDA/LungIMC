library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
# cell_spatial_path <- file.path(phenotype_dir, paste0("InteractionsTestIter200", ".RData"))
cell_spatial_path <- file.path(phenotype_dir, "DelaunayInteraction", paste0("DelaunayInteractionThreshold50", ".RData"))
load(cell_spatial_path)

# replace NA to -1, both unavailable 0
cell_type_num <- 15
interaction_num <- cell_type_num * cell_type_num
num_roi <- interaction_out@nrows / interaction_num
for (nn in 1:num_roi) {
    start_ind <- (nn - 1) * interaction_num + 1
    end_ind <- nn * interaction_num
    cur_sigvals <- interaction_out$sigval[start_ind:end_ind]
    cell_exists <- rep(0, cell_type_num)
    for (mm in 1:cell_type_num) {
        cur_cell_vals <- cur_sigvals[(mm-1)*cell_type_num+1:mm*cell_type_num]
        if (all(is.na(cur_cell_vals)))
            cell_exists[mm] <- 1
    }
    cur_sigvals <- replace(cur_sigvals, is.na(cur_sigvals), -1)
    missing_inds <- which(cell_exists == 1)
    if (length(missing_inds) > 0) {
        ind_combs <- crossing(missing_inds, missing_inds)
        for (jj in 1:dim(ind_combs)[1]) {
            ele_ind <- (ind_combs[[jj,1]]-1)* cell_type_num + ind_combs[[jj,2]]
            cur_sigvals[ele_ind] <- 0
        }
    }
    interaction_out$sigval[start_ind:end_ind] <- cur_sigvals
}

## load ROI diagnosis information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_path <- file.path(metadata_dir, "ROI_Info.xlsx")
roi_meta_info <- read.xlsx(roi_info_path)

# AdjacentNormal of Normal/AAH/AIS/MIA/ADC
path_stage <- "ADC"
subset_roi_info <- subset(roi_meta_info, ROI_Diag==path_stage & ROI_Location=="AdjacentNormal")


subset_roi_lst <- subset_roi_info$ROI_ID
subset_roi_num <- length(subset_roi_lst)
subset_out <- interaction_out[interaction_out$group_by %in% subset_roi_lst, ]

from_order <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell", "CD8-T-Cell", "CD4-T-Cell", 
                "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte", "MDSC", "Fibroblast", "Undefined")
to_order <- c("Undefined", "Fibroblast", "MDSC", "Monocyte", "Macrophage", "Proliferating-Cell", "T-Reg-Cell", "CD4-T-Cell", "CD8-T-Cell", 
              "Endothelial-Cell", "Dendritic-Cell", "NK-Cell", "Neutrophil", "B-Cell", "Epithelial-Cell")

max_per_val <- 1.000
min_per_val <- -1.000

subset_out %>% as_tibble() %>%
    group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / subset_roi_num) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val))) %>%
    mutate(from_label=factor(from_label, levels=from_order)) %>%
    mutate(to_label=factor(to_label, levels=to_order)) %>%
    ggplot() +
    geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
    scale_fill_continuous(limits=c(min_per_val, max_per_val), breaks=seq(5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))