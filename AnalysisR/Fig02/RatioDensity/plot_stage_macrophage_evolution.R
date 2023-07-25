library(tidyverse)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
subtype_proportion_dir <- file.path(data_root_dir, "NatureFigures", "Fig02", "ProportionDensity", "ProportionEvolution")
stage_subcell_ratio_path <- file.path(subtype_proportion_dir, "stage_macrophage_ratios.csv")
cell_ratio_df <- read.csv(stage_subcell_ratio_path)

cell_ratio_df$Stage <- factor(cell_ratio_df$Stage, 
                              levels = c("Normal", "AAH", "AIS", "MIA", "ADC"))
cell_type_orders <- c("CD163- Macrophage", "CD163- PDL1+ Macrophage", "CD163+ Ki67+ Macrophage", "CD163+ Macrophage")
cell_ratio_df$CellType <- factor(cell_ratio_df$CellType, levels = rev(cell_type_orders))
ggplot(cell_ratio_df, aes(x = Stage, y = CellRatio, fill = CellType)) + 
    geom_bar(stat="identity", position = "fill") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Cell Subtype Cell Proportions")

stage_all_cell_plot_path <- file.path(subtype_proportion_dir, "stage_macrophage_ratios.pdf")
ggsave(filename = stage_all_cell_plot_path, device='pdf', width=10, height=8, dpi=300)
while (!is.null(dev.list()))
    dev.off()    