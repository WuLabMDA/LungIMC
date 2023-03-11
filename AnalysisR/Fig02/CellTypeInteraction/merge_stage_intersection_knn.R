library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
cell_spatial_path <- file.path(phenotype_dir, paste0("InteractionsTestIter200KNN", ".RData"))
load(cell_spatial_path)

## load ROI diagnosis information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_path <- file.path(metadata_dir, "ROI_Info.xlsx")
roi_meta_info <- read.xlsx(roi_info_path)


from_order <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell", "CD8-T-Cell", "CD4-T-Cell", 
                "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte", "MDSC", "Fibroblast", "Undefined")
to_order <- c("Undefined", "Fibroblast", "MDSC", "Monocyte", "Macrophage", "Proliferating-Cell", "T-Reg-Cell", "CD4-T-Cell", "CD8-T-Cell", 
              "Endothelial-Cell", "Dendritic-Cell", "NK-Cell", "Neutrophil", "B-Cell", "Epithelial-Cell")


# boundary-4
max_per_val <- 1.0
min_per_val <- -1.0

# subset roi information
normal_roi_info <- subset(roi_meta_info, ROI_Diag=="Normal")
aah_roi_info <- subset(roi_meta_info, ROI_Diag=="AAH" & ROI_Location=="Tumor")
ais_roi_info <- subset(roi_meta_info, ROI_Diag=="AIS" & ROI_Location=="Tumor")
mia_roi_info <- subset(roi_meta_info, ROI_Diag=="MIA" & ROI_Location=="Tumor")
adc_roi_info <- subset(roi_meta_info, ROI_Diag=="ADC" & ROI_Location=="Tumor")

# subset interaction
normal_subset_out <- interaction_out[interaction_out$group_by %in% normal_roi_info$ROI_ID, ]
aah_subset_out <- interaction_out[interaction_out$group_by %in% aah_roi_info$ROI_ID, ]
ais_subset_out <- interaction_out[interaction_out$group_by %in% ais_roi_info$ROI_ID, ]
mia_subset_out <- interaction_out[interaction_out$group_by %in% mia_roi_info$ROI_ID, ]
adc_subset_out <- interaction_out[interaction_out$group_by %in% adc_roi_info$ROI_ID, ]


# update normal
normal_subset <- normal_subset_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(normal_roi_info$ROI_ID)) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
normal_subset$to_label <- paste0(normal_subset$to_label, "-Normal")
# update aah
aah_subset <- aah_subset_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(aah_roi_info$ROI_ID)) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
aah_subset$to_label <- paste0(aah_subset$to_label, "-AAH")
# update ais
ais_subset <- ais_subset_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(ais_roi_info$ROI_ID)) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
ais_subset$to_label <- paste0(ais_subset$to_label, "-AIS")
# update mia
mia_subset <- mia_subset_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(mia_roi_info$ROI_ID)) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
mia_subset$to_label <- paste0(mia_subset$to_label, "-MIA")
# update mia
adc_subset <- adc_subset_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(adc_roi_info$ROI_ID)) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
adc_subset$to_label <- paste0(adc_subset$to_label, "-ADC")



merge_subset <- rbind(normal_subset, aah_subset, ais_subset, mia_subset, adc_subset)
merge_to_order <- c()
for (cell_type in to_order) 
    for (stage_type in c("ADC", "MIA", "AIS", "AAH", "Normal")) 
        merge_to_order <- append(merge_to_order, paste(cell_type, stage_type, sep="-"))


order_merge_set <- merge_subset %>% as_tibble() %>% 
    mutate(from_label=factor(from_label, levels=from_order)) %>%
    mutate(to_label=factor(to_label, levels=merge_to_order))

order_merge_set %>%
    ggplot() +
    geom_tile(aes(from_label, to_label, fill=sum_sigval)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
    theme(axis.text.x=element_text(angle=45, hjust=1)) 
