library(tidyverse)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
subtype_proportion_dir <- file.path(data_root_dir, "NatureFigures", "Fig02", "ProportionDensity", "ProportionEvolution")
stage_subcell_ratio_path <- file.path(subtype_proportion_dir, "stage_subtype_cell_ratios.csv")
cell_ratio_df <- read.csv(stage_subcell_ratio_path)

cell_ratio_df$Stage <- factor(cell_ratio_df$Stage, 
                              levels = c("Normal", "AAH", "AIS", "MIA", "ADC"))
cell_type_orders <- c("Ki67+ Epithelial", "Ki67+ PDL1+ Epithelial", "CD73+ Epithelial", "Other Epithelial",
                      "Ki67+ B-Cells", "Ki67- B-Cells", "Neutrophil", "Ki67+ NK Cells", "Ki67- NK Cells", 
                      "Ki67+ Dendritic Cells", "HLADR+ Dendritic Cells", "Other Dendritic Cells", "Endothelial-Cell",
                      "Cytotoxic CD8 T-Cells", "Memory CD8 T-Cells", "Exhausted CD8 T-Cells", "Ki67+ CD8 T-Cells", "Naive CD8 T-Cells",
                      "Memory CD4 T-Cells", "Exhausted CD4 T-Cells", "Ki67+ CD4 T-Cells", "Naive CD4 T-Cells", 
                      "Ki67+ Treg-Cells", "Ki67- Treg-Cells", "Proliferating-Cell", "CD163+ Macrophage",
                      "CD163+ Ki67+ Macrophage", "CD163- Macrophage", "CD163- PDL1+ Macrophage", "Monocyte",
                      "MDSC", "Fibroblast", "Undefined")
cell_ratio_df$CellType <- factor(cell_ratio_df$CellType, levels = rev(cell_type_orders))
ggplot(cell_ratio_df, aes(x = Stage, y = CellRatio, fill = CellType)) + 
    geom_bar(stat="identity", position = "fill") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Cell Subtype Cell Proportions")

stage_all_cell_plot_path <- file.path(subtype_proportion_dir, "stage_subtype_cell_ratios.pdf")
ggsave(filename = stage_all_cell_plot_path, device='pdf', width=10, height=8, dpi=300)
while (!is.null(dev.list()))
    dev.off()    