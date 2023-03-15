library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
cell_spatial_path <- file.path(phenotype_dir, paste0("InteractionsTestIter200Delaunay", ".RData"))
load(cell_spatial_path)

## load ROI diagnosis information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_path <- file.path(metadata_dir, "ROI_Info.xlsx")
roi_meta_info <- read.xlsx(roi_info_path)

## load tbm information
slide_info_path <- file.path(metadata_dir, "TMB", "LungSlideTMB2.csv")
roi_slide_info <- read.csv(slide_info_path)

# AAH/AIS/MIA/ADC
path_stage <- "ADC"
subset_roi_info <- subset(roi_meta_info, ROI_Diag==path_stage & ROI_Location=="Tumor")
roi_slides <- str_extract_all(subset_roi_info$ROI_ID, ".+(?=-ROI)", simplify = TRUE)
subset_roi_info <- cbind(subset_roi_info, roi_slides)

subset_smoke_info <- filter(roi_slide_info, Stages==path_stage)
low_subset_smoke <- filter(subset_smoke_info, Smoke=="Non-smoker")
high_subset_smoke <- filter(subset_smoke_info, Smoke=="Smoker")


low_sub_roi_info <- filter(subset_roi_info, roi_slides %in% low_subset_smoke$Slides)
high_sub_roi_info <- filter(subset_roi_info, roi_slides %in% high_subset_smoke$Slides)
low_sub_roi_lst <- low_sub_roi_info$ROI_ID
high_sub_roi_lst <- high_sub_roi_info$ROI_ID

from_order <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell", "CD8-T-Cell", "CD4-T-Cell", 
                "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte", "MDSC", "Fibroblast", "Undefined")
to_order <- c("Undefined", "Fibroblast", "MDSC", "Monocyte", "Macrophage", "Proliferating-Cell", "T-Reg-Cell", "CD4-T-Cell", "CD8-T-Cell", 
              "Endothelial-Cell", "Dendritic-Cell", "NK-Cell", "Neutrophil", "B-Cell", "Epithelial-Cell")

low_subset_out <- interaction_out[interaction_out$group_by %in% low_sub_roi_lst, ]
high_subset_out <- interaction_out[interaction_out$group_by %in% high_sub_roi_lst, ]

max_per_val <- 1.00
min_per_val <- -0.6


low_subset <- low_subset_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(low_sub_roi_lst)) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
low_subset$from_label <- paste0(low_subset$from_label, "-Non-smoker")

high_subset <- high_subset_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(high_sub_roi_lst)) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
high_subset$from_label <- paste0(high_subset$from_label, "-Smoker")

merge_subset <- rbind(low_subset, high_subset)

merge_from_order <- c()
for (cell_type in from_order) 
    for (tmb_type in c("Non-smoker", "Smoker")) 
        merge_from_order <- append(merge_from_order, paste(cell_type, tmb_type, sep="-"))


order_merge_set <- merge_subset %>% as_tibble() %>% 
    mutate(from_label=factor(from_label, levels=merge_from_order)) %>%
    mutate(to_label=factor(to_label, levels=to_order))


order_merge_set %>%
    ggplot() +
    geom_tile(aes(from_label, to_label, fill=sum_sigval)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
    theme(axis.text.x=element_text(angle=45, hjust=1)) 

