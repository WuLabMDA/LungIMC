library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)
library(NCmisc)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
# cell_spatial_path <- file.path(phenotype_dir, paste0("InteractionsTestIter200", ".RData"))
cell_spatial_path <- file.path(phenotype_dir, "DelaunayInteraction", paste0("DelaunayInteractionThreshold50", ".RData"))
load(cell_spatial_path)

## load ROI diagnosis information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_path <- file.path(metadata_dir, "ROI_Info.xlsx")
roi_meta_info <- read.xlsx(roi_info_path)

lesion_info_path <- file.path(metadata_dir, "Lesion_Info.xlsx")
lesion_meta_info <- read.xlsx(lesion_info_path)
patient_info_path <- file.path(metadata_dir, "Patient_Info.xlsx")
patient_meta_info <- read.xlsx(patient_info_path)
lesion_smoke_df <- left_join(lesion_meta_info, patient_meta_info, by = join_by(Patient_ID == PatientID))
lesion_normal_df <- subset(lesion_smoke_df, Slide_Diag == "Normal")
low_subset_smoke <- subset(lesion_normal_df, SmokeStatus=="Non-Smoker")
high_subset_smoke <- subset(lesion_normal_df, SmokeStatus=="Smoker")
subset_roi_info <- subset(roi_meta_info, ROI_Diag=="Normal")
roi_slides <- str_extract_all(subset_roi_info$ROI_ID, ".+(?=-ROI)", simplify = TRUE)
subset_roi_info <- cbind(subset_roi_info, roi_slides)
low_sub_roi_info <- filter(subset_roi_info, roi_slides %in% low_subset_smoke$Slide_ID)
high_sub_roi_info <- filter(subset_roi_info, roi_slides %in% high_subset_smoke$Slide_ID)
low_sub_roi_lst <- low_sub_roi_info$ROI_ID
high_sub_roi_lst <- high_sub_roi_info$ROI_ID


from_order <- c("Epithelial-Cell", "B-Cell", "Neutrophil", "NK-Cell", "Dendritic-Cell", "Endothelial-Cell", "CD8-T-Cell", "CD4-T-Cell", 
                "T-Reg-Cell", "Proliferating-Cell", "Macrophage", "Monocyte", "MDSC", "Fibroblast", "Undefined")
to_order <- c("Undefined", "Fibroblast", "MDSC", "Monocyte", "Macrophage", "Proliferating-Cell", "T-Reg-Cell", "CD4-T-Cell", "CD8-T-Cell", 
              "Endothelial-Cell", "Dendritic-Cell", "NK-Cell", "Neutrophil", "B-Cell", "Epithelial-Cell")

low_subset_out <- interaction_out[interaction_out$group_by %in% low_sub_roi_lst, ]
high_subset_out <- interaction_out[interaction_out$group_by %in% high_sub_roi_lst, ]


low_subset <- low_subset_out %>% as_tibble() %>% group_by(from_label, to_label)
high_subset <- high_subset_out %>% as_tibble() %>% group_by(from_label, to_label)


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
        low_pair <- low_subset[low_subset$from_label == p_from_label & low_subset$to_label == p_to_label,] %>% drop_na(sigval)
        high_pair <- high_subset[high_subset$from_label == p_from_label & high_subset$to_label == p_to_label,] %>% drop_na(sigval)
        if (length(low_pair$sigval) <= 3 | length(high_pair$sigval) <= 3)
            cur_p_val <- 0.99
        else {
            merge_sig_vals <- sigvals <- c(high_pair$sigval, low_pair$sigval)
            if (all(merge_sig_vals == merge_sig_vals[1]))
                cur_p_val <- 0.99
            else {
                test_pval <- t.test(high_pair$sigval, low_pair$sigval)
                cur_p_val <- test_pval$p.value
            }            
        }
        pval_df[p_to_label, p_from_label] <- cur_p_val
        group_pvals[group_ind, ] <- c(p_to_label, p_from_label, cur_p_val)
        group_ind <- group_ind + 1
    }
}


adjust_group <- group_pvals %>% 
    mutate(pVal = as.numeric(pVal)) %>%
    mutate(gGroup = case_when(pVal > 0.05 ~ 'NS', pVal > 0.01 ~ '*', pVal > 0.001 ~ '**', .default = "***"))


adjust_group %>% ggplot() +
    geom_point(aes(x = factor(ToType, level=from_order), y = factor(FromType, level=to_order), 
                   size = p.to.Z(pVal), col = factor(gGroup, level=c("NS", "*", "**", "***")))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

