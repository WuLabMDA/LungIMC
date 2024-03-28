library(tidyverse)

# prepare directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
tim3_dir <- file.path(phenotype_dir, "TIM3")
tim3_path <- file.path(tim3_dir, "cell_tim3.rds")
cell_exp_df <- readRDS(tim3_path)
tim3_roi_proportion_dir <- file.path(tim3_dir, "CellType-ROI")
if (!file.exists(tim3_roi_proportion_dir))
  dir.create(tim3_roi_proportion_dir, recursive = TRUE)

# calculate stage-wise mean expression
cell_stages <- cell_exp_df$cell_stage
cell_PDL1_vals <- cell_exp_df$PDL1

PDL1_stage_exprs <- c()
unique_stages <- c("Normal", "AAH", "AIS", "MIA", "ADC")
for (cur_stage in unique_stages) {
  cur_stage_indices <- which(cell_stages == cur_stage)
  cur_PDL1_vals <- cell_PDL1_vals[cur_stage_indices]
  PDL1_stage_exprs[cur_stage] <- sum(cur_PDL1_vals) / length(cur_PDL1_vals)
}

# plot
barplot(PDL1_stage_exprs)