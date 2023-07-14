library(tidyverse)
library(ggpubr)

# load data
data_root_dir <- "E:/LungIMCData/HumanWholeIMC/Aggregation/MarginCore"
adc_loc_fea_path <- file.path(data_root_dir, "adc_lesion_info.csv")
adc_loc_fea <- read.csv(adc_loc_fea_path)

# feature list
study_fea_lst = c("CD8-T-Cell-Proportion", "CD8-T-Cell-Density", "Macrophage-Proportion", "Macrophage-Density", 
                 "B-Cell-Proportion", "B-Cell-Density", "Endothelial-Cell-Proportion", "Endothelial-Cell-Density",
                 "CD8-T-Cell-Epithelial-Cell", "Macrophage-Epithelial-Cell", "Endothelial-Cell-Fibroblast", "Neutrophil-T-Reg-Cell")


for (cur_fea in study_fea_lst) {
    fea_name <- gsub("-", ".", cur_fea)
    margin_fea_df <- filter(adc_loc_fea, Loc == "TumorMargin")
    core_fea_df <- filter(adc_loc_fea, Loc == "TumorCore")
    margin_fea <- margin_fea_df[, fea_name]
    core_fea <- core_fea_df[, fea_name]
    cmp_df <- data.frame(Margin=margin_fea, Core=core_fea)
    
    # plot
    p <- ggpaired(cmp_df, cond1 = "Margin", cond2 = "Core", fill = "condition", 
             line.color = "gray", line.size = 0.4, palette = "npg") + ylab(fea_name)
    cmp_plot_path <- file.path(data_root_dir, paste(cur_fea, "png", sep = "."))
    png(cmp_plot_path)
    print(p)
    dev.off()
}




