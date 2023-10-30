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

# extract cell_types & TIM3 expression
cell_types <- cell_exp_df$cell_type
cell_rois <- cell_exp_df$cell_roi
cell_stages <- cell_exp_df$cell_stage
cell_tim3_vals <- cell_exp_df$TIM3

unique_cell_types <- c("Epithelial-Cell", "Endothelial-Cell", "Fibroblast", "CD4-T-Cell", "CD8-T-Cell", 
                       "T-Reg-Cell", "B-Cell", "Macrophage", "Monocyte", "Dendritic-Cell", 
                       "Neutrophil", "MDSC", "NK-Cell", "Proliferating-Cell")
for (cur_type in unique_cell_types) {
    cur_cell_indices <- which(cell_types == cur_type)
    cur_cell_rois <- cell_rois[cur_cell_indices]
    cur_cell_stages <- cell_stages[cur_cell_indices]
    cur_tim3_vals <- cell_tim3_vals[cur_cell_indices]
    
    # find unique rois
    unique_rois <- unique(cur_cell_rois)
    stage_lst <- c()
    tim3_ratio_lst <- c()
    for (cur_roi in unique_rois) {
        cur_roi_indices <- which(cur_cell_rois == cur_roi)
        cur_roi_stages <- cur_cell_stages[cur_roi_indices]
        if (length(unique(cur_roi_stages)) != 1)
            print("Multiple stages")
        stage_lst <- append(stage_lst, cur_roi_stages[1])
        cur_roi_tim3_vals <- cur_tim3_vals[cur_roi_indices]
        cur_roi_tim3_postive_ratio <- sum(cur_roi_tim3_vals > 0.0) / length(cur_roi_tim3_vals)
        tim3_ratio_lst <- append(tim3_ratio_lst, cur_roi_tim3_postive_ratio)
    }
    roi_ratio_df <- data.frame(Proportion=tim3_ratio_lst, Stage=stage_lst)
    roi_ratio_df$Stage <- factor(roi_ratio_df$Stage, levels = c("Normal", "AAH", "AIS", "MIA", "ADC"))
    ggplot(roi_ratio_df, aes(x = Stage, y = Proportion, colour = Stage)) + 
        geom_boxplot() +
        geom_jitter(width = 0.1) + 
        ylim(0.0, 1.0) +
        theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        ggtitle(paste0(cur_type, " - Proportion"))
    
    tim3_proportion_plot_path <- file.path(tim3_roi_proportion_dir, paste0(cur_type, ".pdf"))
    ggsave(filename = tim3_proportion_plot_path, device='pdf', width=5, height=8, dpi=300)
    while (!is.null(dev.list()))
        dev.off()
}