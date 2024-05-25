library(tidyverse)

# prepare directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
# data_root_dir <- "/Volumes/wulab/Ping/HumanIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
tim3_dir <- file.path(phenotype_dir, "TIM3")
tim3_path <- file.path(tim3_dir, "all_cell_tim3.rds")
cell_exp_df <- readRDS(tim3_path)
cell_tim3_vals <- cell_exp_df$TIM3

# calculate stage-wise mean expression
cell_locations <- cell_exp_df$cell_location
cell_stages <- cell_exp_df$cell_stage

for (cell_ind in 1:length(cell_exp_df$cell_id)) {
  if (cell_locations[cell_ind] == "Normal")
    cell_stages[cell_ind] <- "Normal"
}
cell_exp_df$cell_stage <- cell_stages

# # violin plot 
# cell_exp_df$cell_stage <- factor(cell_exp_df$cell_stage, levels=c("Normal", "AAH", "AIS", "MIA", "ADC"))
# ggplot(cell_exp_df, aes(x = cell_stage, y = TIM3)) + geom_violin() + ylim(0.01, 0.16) + geom_jitter(size = 0.2, width = 0.3)
# tim3_expression_plot_path <- file.path(tim3_dir, "TIMC3-Expression-Violin-v1.pdf")
# ggsave(filename = tim3_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
# while (!is.null(dev.list()))
#   dev.off()

unique_stages <- c("Normal", "AAH", "AIS", "MIA", "ADC")
for (cur_stage in unique_stages) {
    cur_stage_indices <- which(cell_stages == cur_stage)
    cur_tim3_vals <- cell_tim3_vals[cur_stage_indices]
    tim3_stage_exprs[cur_stage] <- sum(cur_tim3_vals) / length(cur_tim3_vals)
}
# barplot
barplot(tim3_stage_exprs)
