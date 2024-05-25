library(tidyverse)

# prepare directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
tim3_dir <- file.path(phenotype_dir, "TIM3")
tim3_path <- file.path(tim3_dir, "all_cell_tim3.rds")
cell_exp_df <- readRDS(tim3_path)

# extract cell_types & TIM3 expression
cell_types <- cell_exp_df$cell_type
cell_tim3_vals <- cell_exp_df$TIM3

unique_cell_types <- c("Epithelial-Cell", "Endothelial-Cell", "Fibroblast", "CD4-T-Cell", "CD8-T-Cell", 
                       "T-Reg-Cell", "B-Cell", "Macrophage", "Monocyte", "Dendritic-Cell", 
                       "Neutrophil", "MDSC", "NK-Cell", "Proliferating-Cell")

# calculate TIM3 postivity ratios
tim3_positive_ratios <- c()
for (ind in 1:length(unique_cell_types)) {
    cur_cell_type <- unique_cell_types[ind]
    cur_cell_indices <- which(cur_cell_type == cell_types)
    cur_tim3_vals <- cell_tim3_vals[cur_cell_indices]
    cur_tim3_postive_num <- sum(cur_tim3_vals > 0.0)
    cur_tim3_postive_ratio <- cur_tim3_postive_num * 1.0 / length(cur_cell_indices)
    tim3_positive_ratios[cur_cell_type] <- cur_tim3_postive_ratio
}

# plot
barplot(tim3_positive_ratios, main = "Major Cell Type TIM3 Positivity Proportion", 
        xlab = "Cell Types", ylab = "Proportion", cex.names=0.5)