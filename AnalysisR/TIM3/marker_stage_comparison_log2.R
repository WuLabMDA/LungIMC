library(tidyverse)

# prepare directory
# data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
data_root_dir <- "/Volumes/wulab/Ping/HumanIMCData/HumanWholeIMC"

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

# ROI levels
stage_lst <-c()
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
  
  for (cur_roi in unique_rois) {
    cur_cell_indices <- which(cur_stage_rois == cur_roi)
    stage_lst <- append(stage_lst, cur_stage)
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
}

roi_marker_exp_df <- data.frame(Stage = stage_lst,
                                TIM3 = tim3_lst, ICOS = icos_lst, TIGIT = tigit_lst, CD73 = cd73_lst,
                                PDL1 = pdl1_lst, LAG3 = lag3_lst, VISTA = vista_lst, PD1 = pd1_lst,
                                B7H3 = b7h3_lst, CTLA4 = ctla4_lst, IDO1 = ido1_lst)
roi_marker_exp_df$Stage <- factor(roi_marker_exp_df$Stage, levels=c("Normal", "AAH", "AIS", "MIA", "ADC"))

# TIM3
ggplot(roi_marker_exp_df, aes(x = Stage, y = log2(TIM3))) + geom_boxplot() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)
# add pairwise.wilcox.test test
# test_res <- pairwise.wilcox.test(roi_marker_exp_df$TIM3, roi_marker_exp_df$Stage, p.adjust.method = "bonferroni") 

pval <- kruskal.test(TIM3 ~ Stage, data = roi_marker_exp_df) 
marker_expression_plot_path <- file.path(marker_expression_cmp_dir, paste0("TIM3-Stage-BoxplotLog2-pval-", toString(pval$p.value), ".pdf"))
ggsave(filename = marker_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
while (!is.null(dev.list()))
  dev.off()

# ICOS
ggplot(roi_marker_exp_df, aes(x = Stage, y = log2(ICOS))) + geom_boxplot() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)
pval <- kruskal.test(ICOS ~ Stage, data = roi_marker_exp_df) 
marker_expression_plot_path <- file.path(marker_expression_cmp_dir, paste0("ICOS-Stage-BoxplotLog2-pval-", toString(pval$p.value), ".pdf"))
ggsave(filename = marker_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
while (!is.null(dev.list()))
  dev.off()

# TIGIT
ggplot(roi_marker_exp_df, aes(x = Stage, y = log2(TIGIT))) + geom_boxplot() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)
pval <- kruskal.test(TIGIT ~ Stage, data = roi_marker_exp_df) 
marker_expression_plot_path <- file.path(marker_expression_cmp_dir, paste0("TIGIT-Stage-BoxplotLog2-pval-", toString(pval$p.value), ".pdf"))
ggsave(filename = marker_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
while (!is.null(dev.list()))
  dev.off()

# CD73
ggplot(roi_marker_exp_df, aes(x = Stage, y = log2(CD73))) + geom_boxplot() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)
pval <- kruskal.test(CD73 ~ Stage, data = roi_marker_exp_df) 
marker_expression_plot_path <- file.path(marker_expression_cmp_dir, paste0("CD73-Stage-BoxplotLog2-pval-", toString(pval$p.value), ".pdf"))
ggsave(filename = marker_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
while (!is.null(dev.list()))
  dev.off()

# PDL1
ggplot(roi_marker_exp_df, aes(x = Stage, y = log2(PDL1))) + geom_boxplot() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)
pval <- kruskal.test(PDL1 ~ Stage, data = roi_marker_exp_df) 
marker_expression_plot_path <- file.path(marker_expression_cmp_dir, paste0("PDL1-Stage-BoxplotLog2-pval-", toString(pval$p.value), ".pdf"))
ggsave(filename = marker_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
while (!is.null(dev.list()))
  dev.off()

# LAG3
ggplot(roi_marker_exp_df, aes(x = Stage, y = log2(LAG3))) + geom_boxplot() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)
pval <- kruskal.test(LAG3 ~ Stage, data = roi_marker_exp_df) 
marker_expression_plot_path <- file.path(marker_expression_cmp_dir, paste0("LAG3-Stage-BoxplotLog2-pval-", toString(pval$p.value), ".pdf"))
ggsave(filename = marker_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
while (!is.null(dev.list()))
  dev.off()

# VISTA
ggplot(roi_marker_exp_df, aes(x = Stage, y = log2(VISTA))) + geom_boxplot() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)
pval <- kruskal.test(VISTA ~ Stage, data = roi_marker_exp_df) 
marker_expression_plot_path <- file.path(marker_expression_cmp_dir, paste0("VISTA-Stage-BoxplotLog2-pval-", toString(pval$p.value), ".pdf"))
ggsave(filename = marker_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
while (!is.null(dev.list()))
  dev.off()

# PD1
ggplot(roi_marker_exp_df, aes(x = Stage, y = log2(PD1))) + geom_boxplot() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)
pval <- kruskal.test(PD1 ~ Stage, data = roi_marker_exp_df) 
marker_expression_plot_path <- file.path(marker_expression_cmp_dir, paste0("PD1-Stage-BoxplotLog2-pval-", toString(pval$p.value), ".pdf"))
ggsave(filename = marker_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
while (!is.null(dev.list()))
  dev.off()

# B7H3
ggplot(roi_marker_exp_df, aes(x = Stage, y = log2(B7H3))) + geom_boxplot() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)
pval <- kruskal.test(B7H3 ~ Stage, data = roi_marker_exp_df) 
marker_expression_plot_path <- file.path(marker_expression_cmp_dir, paste0("B7H3-Stage-BoxplotLog2-pval-", toString(pval$p.value), ".pdf"))
ggsave(filename = marker_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
while (!is.null(dev.list()))
  dev.off()

# CTLA4
ggplot(roi_marker_exp_df, aes(x = Stage, y = log2(CTLA4))) + geom_boxplot() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)
pval <- kruskal.test(CTLA4 ~ Stage, data = roi_marker_exp_df) 
marker_expression_plot_path <- file.path(marker_expression_cmp_dir, paste0("CTLA4-Stage-BoxplotLog2-pval-", toString(pval$p.value), ".pdf"))
ggsave(filename = marker_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
while (!is.null(dev.list()))
  dev.off()

# IDO1
ggplot(roi_marker_exp_df, aes(x = Stage, y = log2(IDO1))) + geom_boxplot() + ylim(-6.0, 1.0) + geom_jitter(size = 0.2, width = 0.3)
pval <- kruskal.test(IDO1 ~ Stage, data = roi_marker_exp_df) 
marker_expression_plot_path <- file.path(marker_expression_cmp_dir, paste0("IDO1-Stage-BoxplotLog2-pval-", toString(pval$p.value), ".pdf"))
ggsave(filename = marker_expression_plot_path, device='pdf', width=8, height=5, dpi=300)
while (!is.null(dev.list()))
  dev.off()
