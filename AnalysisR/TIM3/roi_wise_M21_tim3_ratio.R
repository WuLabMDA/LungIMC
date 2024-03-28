library(readxl)
library(stringr)
library(tidyverse)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
tim3_dir <- file.path(phenotype_dir, "TIM3")
tim3_path <- file.path(tim3_dir, "cell_subtype_tim3.rds")

cell_exp_df <- readRDS(tim3_path)

# extract cell_types & TIM3 expression
cell_types <- cell_exp_df$cell_type
cell_stages <- cell_exp_df$cell_stage
cell_tim3_vals <- cell_exp_df$TIM3

# filtering positive tim3
tim3_pos_indices <- which(cell_tim3_vals > 0)
tim3_pos_cell_types <- cell_types[tim3_pos_indices]
tim3_pos_cell_stage <- cell_stages[tim3_pos_indices]

# m2/m1 cell types
m2_cell_types <- c("CD163+ Macrophages", "Ki67+ Macrophages")
m1_cell_types <- c("CD163- Macrophages")

M21_ratios <- c()
unique_stages <- c("Normal", "AAH", "AIS", "MIA", "ADC")
for (cur_stage in unique_stages) {
  cur_stage_indices <- which(tim3_pos_cell_stage == cur_stage)
  cur_cell_types <- tim3_pos_cell_types[cur_stage_indices]
  m2_cell_num <- sum(cur_cell_types %in% m2_cell_types)
  m1_cell_num <- sum(cur_cell_types %in% m1_cell_types)
  M21_ratios[cur_stage] <- (m2_cell_num + 0.01) / (m1_cell_num + 0.01) 
  print(paste0(cur_stage, " M2:", m2_cell_num, " M1:", m1_cell_num))
}

# plot
barplot(M21_ratios)