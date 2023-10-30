library(tidyverse)

# prepare directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
tim3_dir <- file.path(phenotype_dir, "TIM3")
tim3_path <- file.path(tim3_dir, "cell_tim3.rds")
cell_exp_df <- readRDS(tim3_path)

# calculate stage-wise mean expression
cell_stages <- cell_exp_df$cell_stage
cell_tim3_vals <- cell_exp_df$TIM3

tim3_stage_exprs <- c()
unique_stages <- c("Normal", "AAH", "AIS", "MIA", "ADC")
for (cur_stage in unique_stages) {
    cur_stage_indices <- which(cell_stages == cur_stage)
    cur_tim3_vals <- cell_tim3_vals[cur_stage_indices]
    tim3_stage_exprs[cur_stage] <- sum(cur_tim3_vals) / length(cur_tim3_vals)
}

# plot
barplot(tim3_stage_exprs)
