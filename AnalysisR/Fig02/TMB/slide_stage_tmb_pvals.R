library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
celltype_delaunay_dir <- file.path(phenotype_dir, "DelaunayInteraction")
threshold_val <- 50
cell_type_interaction_name <- paste0("DelaunayInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_delaunay_dir, cell_type_interaction_name)
load(cell_type_interaction_path)

# load ROI information
metadata_dir <- file.path(data_root_dir, "Metadata")
slide_info_name <- "Lesion_Info_Aggregation.csv"
slide_info_path <- file.path(metadata_dir, slide_info_name)
slide_info_df <- read_csv(slide_info_path)

# AAH/AIS/MIA/ADC
path_stage <- "ADC"
subset_slide_df <- slide_info_df[slide_info_df$Slide_Diag==path_stage,]
subset_slide_df <- subset_slide_df[!is.na(subset_slide_df$TMB),]
subset_slide_lst <- subset_slide_df$Slide_ID
low_slide_df <- subset(subset_slide_df, TMB == "Low")
high_slide_df <- subset(subset_slide_df, TMB == "High")
low_slide_lst <- low_slide_df$Slide_ID
high_slide_lst <- high_slide_df$Slide_ID

## find subset ROIs
roi_info_path <- file.path(metadata_dir, "ROI_Info.xlsx")
roi_meta_info <- read.xlsx(roi_info_path)
subset_roi_info <- subset(roi_meta_info, ROI_Diag==path_stage & ROI_Location=="Tumor")
subset_roi_slides <- str_extract_all(subset_roi_info$ROI_ID, ".+(?=-ROI)", simplify = TRUE)
subset_roi_slides <- subset_roi_slides[, 1]
subset_roi_info$ROI_Slide <- subset_roi_slides
subset_roi_info <- subset(subset_roi_info, ROI_Slide %in% subset_slide_lst) 

from_order <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell", "CD8-T-Cell", "CD4-T-Cell", 
                "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte", "MDSC", "Fibroblast", "Undefined")
to_order <- c("Undefined", "Fibroblast", "MDSC", "Monocyte", "Macrophage", "Proliferating-Cell", "T-Reg-Cell", "CD4-T-Cell", "CD8-T-Cell", 
              "Endothelial-Cell", "Dendritic-Cell", "NK-Cell", "Neutrophil", "B-Cell", "Epithelial-Cell")

max_per_val <- 1.00
min_per_val <- -0.6

low_tmb_interactions <- matrix(0, nrow=length(from_order)*length(to_order), ncol = length(low_slide_lst))
for (ind in 1:length(low_slide_lst)) {
    cur_slide <- low_slide_lst[ind]
    cur_roi_lst <- subset_roi_info[subset_roi_info$ROI_Slide == cur_slide, "ROI_ID"]
    cur_slide_out <- interaction_out[interaction_out$group_by %in% cur_roi_lst, ]
    cur_slide_interaction <- cur_slide_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
        summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(cur_roi_lst)) %>%
        mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
    low_tmb_interactions[, ind] <- cur_slide_interaction$sum_sigval
}

high_tmb_interactions <- matrix(0, nrow=length(from_order)*length(to_order), ncol = length(high_slide_lst))
for (ind in 1:length(high_slide_lst)) {
    cur_slide <- high_slide_lst[ind]
    cur_roi_lst <- subset_roi_info[subset_roi_info$ROI_Slide == cur_slide, "ROI_ID"]
    cur_slide_out <- interaction_out[interaction_out$group_by %in% cur_roi_lst, ]
    cur_slide_interaction <- cur_slide_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
        summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(cur_roi_lst)) %>%
        mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
    high_tmb_interactions[, ind] <- cur_slide_interaction$sum_sigval
}

p_from_order <- as.character(cur_slide_interaction$from_label)
p_to_order <- as.character(cur_slide_interaction$to_label)

group_pvals <- data.frame(matrix(nrow = length(from_order) * length(from_order), ncol = 3))
colnames(group_pvals) <- c("FromType", "ToType", "pVal")
for (ind in 1:length(p_from_order)) {
    low_slide_vals <- low_tmb_interactions[ind, ]
    high_slide_vals <- high_tmb_interactions[ind, ]
    merge_sig_vals <- c(high_slide_vals, low_slide_vals)
    if (all(merge_sig_vals == merge_sig_vals[1]))
        cur_p_val <- 0.99
    else {
        test_pval <- t.test(high_slide_vals, low_slide_vals)
        cur_p_val <- test_pval$p.value
    }
    group_pvals[ind, ] <- c(p_from_order[ind], p_to_order[ind], cur_p_val)
}

adjust_group <- group_pvals %>% 
    mutate(pVal = as.numeric(pVal)) %>%
    mutate(gGroup = case_when(pVal > 0.05 ~ "NS", pVal > 0.01 ~ "*", .default = "**")) %>%
    mutate(gGroup = factor(gGroup, level=c("NS", "*", "**")))


adjust_group %>% ggplot() +
    geom_point(aes(x = factor(FromType, level=from_order), y = factor(ToType, level=to_order), 
                   size = p.to.Z(pVal), col = gGroup)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))