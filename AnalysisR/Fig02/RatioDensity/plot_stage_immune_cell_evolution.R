library(tidyverse)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
proportion_evolution_dir <- file.path(data_root_dir, "NatureFigures", "Fig02", "ProportionDensity", "ProportionEvolution")
stage_immune_cell_ratio_path <- file.path(proportion_evolution_dir, "stage_immune_cell_ratios.csv")
cell_ratio_df <- read.csv(stage_immune_cell_ratio_path)

cell_ratio_df$Stage <- factor(cell_ratio_df$Stage, 
                              levels = c("Normal", "AAH", "AIS", "MIA", "ADC"))

cell_type_orders <- c("Macrophage", "Neutrophil", "Monocyte", "CD8-T-Cell", "MDSC",  
                      "CD4-T-Cell", "NK-Cell", "Dendritic-Cell", "B-Cell", "T-Reg-Cell")
cell_ratio_df$CellType <- factor(cell_ratio_df$CellType, levels = rev(cell_type_orders))
ggplot(cell_ratio_df, aes(x = Stage, y = CellRatio, fill = CellType)) + 
    geom_bar(stat="identity", position = "fill") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Immune Cell Proportions")

stage_immune_cell_plot_path <- file.path(proportion_evolution_dir, "stage_immune_cell_ratios.pdf")
ggsave(filename = stage_immune_cell_plot_path, device='pdf', width=5, height=8, dpi=300)
while (!is.null(dev.list()))
    dev.off()    