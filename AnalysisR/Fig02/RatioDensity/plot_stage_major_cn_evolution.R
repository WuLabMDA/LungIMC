library(tidyverse)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
proportion_evolution_dir <- file.path(data_root_dir, "NatureFigures", "Fig02", "ProportionDensity", "ProportionEvolution")
stage_all_cell_ratio_path <- file.path(proportion_evolution_dir, "stage_major_cn_ratios.csv")
cell_ratio_df <- read.csv(stage_all_cell_ratio_path)

cell_ratio_df$Stage <- factor(cell_ratio_df$Stage, 
                              levels = c("Normal", "AAH", "AIS", "MIA", "ADC"))
cell_type_orders <- c("Undefined-CN1", "Epithelial1-CN2", "Proliferating-CN3", "Epithelial2-CN4",
                      "Endothelial-CN5", "Fibroblast-CN6", "Macrophage-CN7", "Pan-Immune-CN8")
cell_ratio_df$CellType <- factor(cell_ratio_df$CellType, levels = rev(cell_type_orders))
ggplot(cell_ratio_df, aes(x = Stage, y = CellRatio, fill = CellType)) + 
    geom_bar(stat="identity", position = "fill") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Major Cell Proportions")

stage_all_cell_plot_path <- file.path(proportion_evolution_dir, "stage_major_cn_ratios.pdf")
ggsave(filename = stage_all_cell_plot_path, device='pdf', width=5, height=8, dpi=300)
while (!is.null(dev.list()))
    dev.off()    