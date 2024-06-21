library(tidyverse)

# prepare directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
# data_root_dir <- "/Volumes/wulab/Ping/HumanIMCData/HumanWholeIMC"

phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
tim3_dir <- file.path(phenotype_dir, "TIM3")
tim3_path <- file.path(tim3_dir, "all_cell_tim3.rds")
cell_exp_df <- readRDS(tim3_path)

marker_expression_cmp_dir <- file.path(tim3_dir, "ExpressionCmp")
if (!file.exists(marker_expression_cmp_dir))
    dir.create(marker_expression_cmp_dir, recursive = TRUE)

unique_stages <- c("Normal", "AAH", "AIS", "MIA", "ADC")

cell_locations <- cell_exp_df$cell_location
cell_stages <- cell_exp_df$cell_stage
cell_rois <- cell_exp_df$cell_roi

for (cell_ind in 1:length(cell_exp_df$cell_id)) {
  if (cell_locations[cell_ind] == "Normal")
    cell_stages[cell_ind] <- "Normal"
}


cell_TIM3_vals <- cell_exp_df$TIM3
cell_ICOS_vals <- cell_exp_df$ICOS
cell_TIGIT_vals <- cell_exp_df$TIGIT
cell_CD73_vals <- cell_exp_df$CD73
cell_PDL1_vals <- cell_exp_df$PDL1
cell_LAG3_vals <- cell_exp_df$LAG3
cell_VISTA_vals <- cell_exp_df$VISTA
cell_PD1_vals <- cell_exp_df$PD1
cell_B7H3_vals <- cell_exp_df$B7H3
cell_CTLA4_vals <- cell_exp_df$CTLA4
cell_IDO1_vals <- cell_exp_df$IDO1

for (cur_stage in unique_stages) {
    cur_stage_indices <- which(cell_stages == cur_stage)
    cur_stage_rois <- cell_rois[cur_stage_indices]
    cur_stage_TIM3_vals <- cell_TIM3_vals[cur_stage_indices]
    cur_stage_ICOS_vals <- cell_ICOS_vals[cur_stage_indices]
    cur_stage_TIGIT_vals <- cell_TIGIT_vals[cur_stage_indices]
    cur_stage_CD73_vals <- cell_CD73_vals[cur_stage_indices]
    cur_stage_PDL1_vals <- cell_PDL1_vals[cur_stage_indices]
    cur_stage_LAG3_vals <- cell_LAG3_vals[cur_stage_indices]
    cur_stage_VISTA_vals <- cell_VISTA_vals[cur_stage_indices]
    cur_stage_PD1_vals <- cell_PD1_vals[cur_stage_indices]
    cur_stage_B7H3_vals <- cell_B7H3_vals[cur_stage_indices]
    cur_stage_CTLA4_vals <- cell_CTLA4_vals[cur_stage_indices]
    cur_stage_IDO1_vals <- cell_IDO1_vals[cur_stage_indices]
    unique_rois <- unique(cur_stage_rois)
    
    tim3_lst <- c()
    icos_lst <- c()
    tigit_lst <- c()
    cd73_lst <- c()
    pdl1_lst <- c()
    lag3_lst <- c()
    vista_lst <- c()
    pd1_lst <- c()
    b7h3_lst <- c()
    ctla4_lst <- c()
    ido1_lst <- c()
    for (cur_roi in unique_rois) {
        cur_cell_indices <- which(cur_stage_rois == cur_roi)
        tim3_lst <- append(tim3_lst, sum(cur_stage_TIM3_vals[cur_cell_indices]) / length(cur_cell_indices))
        icos_lst <- append(icos_lst, sum(cur_stage_ICOS_vals[cur_cell_indices]) / length(cur_cell_indices))
        tigit_lst <- append(tigit_lst, sum(cur_stage_TIGIT_vals[cur_cell_indices]) / length(cur_cell_indices))
        cd73_lst <- append(cd73_lst, sum(cur_stage_CD73_vals[cur_cell_indices]) / length(cur_cell_indices))
        pdl1_lst <- append(pdl1_lst, sum(cur_stage_PDL1_vals[cur_cell_indices]) / length(cur_cell_indices))
        lag3_lst <- append(lag3_lst, sum(cur_stage_LAG3_vals[cur_cell_indices]) / length(cur_cell_indices))
        vista_lst <- append(vista_lst, sum(cur_stage_VISTA_vals[cur_cell_indices]) / length(cur_cell_indices))
        pd1_lst <- append(pd1_lst, sum(cur_stage_PD1_vals[cur_cell_indices]) / length(cur_cell_indices))
        b7h3_lst <- append(b7h3_lst, sum(cur_stage_B7H3_vals[cur_cell_indices]) / length(cur_cell_indices))
        ctla4_lst <- append(ctla4_lst, sum(cur_stage_CTLA4_vals[cur_cell_indices]) / length(cur_cell_indices))
        ido1_lst <- append(ido1_lst, sum(cur_stage_IDO1_vals[cur_cell_indices]) / length(cur_cell_indices))
    }
    roi_marker_exp_df <- data.frame(TIM3 = tim3_lst, ICOS = icos_lst, TIGIT = tigit_lst, CD73 = cd73_lst,
                                    PDL1 = pdl1_lst, LAG3 = lag3_lst, VISTA = vista_lst, PD1 = pd1_lst,
                                    B7H3 = b7h3_lst, CTLA4 = ctla4_lst, IDO1 = ido1_lst)
    # ggplot(stack(roi_marker_exp_df), aes(x = ind, y = log10(values))) + geom_violin() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)
    ggplot(stack(roi_marker_exp_df), aes(x = ind, y = log2(values))) + geom_boxplot() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)

    marker_expression_plot_path <- file.path(marker_expression_cmp_dir, paste0(cur_stage, "-BoxLog2.pdf"))
    ggsave(filename = marker_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
    while (!is.null(dev.list()))
        dev.off()
}