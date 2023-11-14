library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
celltype_neighbor_dir <- file.path(phenotype_dir, "DelaunayInteraction")
cell_subtype_interaction_name <- "Subtype23DelaunayInteractionThreshold50.RData"
cell_subtype_interaction_path <- file.path(celltype_neighbor_dir, cell_subtype_interaction_name)
load(cell_subtype_interaction_path)

# replace NA to -1, both unavailable 0
cell_type_num <- 23
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
from_order <- c("Epithelial", "B-Cells", "Neutrophil", "NK-Cells", "Dendritic-Cell", "Endothelial-Cell",
                "Naive CD8 T-Cells", "Cytotoxic CD8 T-Cells", "Memory CD8 T-Cells", "Exhausted CD8 T-Cells", "Ki67+ CD8 T-Cells", 
                "Naive CD4 T-Cells", "Memory CD4 T-Cells", "Exhausted CD4 T-Cells", "Ki67+ CD4 T-Cells", 
                "Treg-Cells", "Proliferating-Cell", "CD163+ Macrophages", "CD163- Macrophages", 
                "Monocyte", "MDSC", "Fibroblast", "Undefined")
to_order <- rev(from_order)
max_per_val <- 1.000
min_per_val <- -1.000

subtype_interaction_dir <- file.path(data_root_dir, "NatureFigures", "Fig02", "AdjacentDistantInteraction", "23AdjacentNormal")
if (!file.exists(subtype_interaction_dir))
    dir.create(subtype_interaction_dir, recursive = TRUE)

# AdjacentNormal of AAH/AIS/MIA/ADC
interested_stages <- c("AAH", "AIS", "MIA", "ADC")
for (path_stage in interested_stages) {
    subset_roi_info <- subset(roi_meta_info, ROI_Diag==path_stage & ROI_Location=="AdjacentNormal")
    subset_roi_lst <- subset_roi_info$ROI_ID
    subset_roi_num <- length(subset_roi_lst)
    subset_out <- interaction_out[interaction_out$group_by %in% subset_roi_lst, ]
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
    subtype_interaction_plot_path <- file.path(subtype_interaction_dir, paste0("Adjacent-", path_stage, ".pdf"))
    ggsave(filename = subtype_interaction_plot_path, device='pdf', width=18, height=15, dpi=300)
    while (!is.null(dev.list()))
        dev.off()    
}