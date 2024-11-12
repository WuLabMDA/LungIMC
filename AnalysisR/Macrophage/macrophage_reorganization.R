library(tidyverse)
library(readxl)
library(stringr)
library(SpatialExperiment)
library(hopkins)
set.seed(1) 

# load data
data_root_dir <- "/Volumes/wulab/Ping/HumanIMCData/HumanWholeIMC"
macrophage_data_path <- file.path(data_root_dir, "Macrophage", "Macrophages.rds")
spe <- readRDS(macrophage_data_path)


# cell id
# cell_ids <- rownames(colData(spe))
cell_rois <- colData(spe)$sample_id
cell_unique_rois <- unique(cell_rois)
# cell_types <- colData(spe)$celltypes
cell_locs <- spatialCoords(spe)

# load ROI information
metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_name <- "ROI_Info_Aggregation.csv"
roi_info_path <- file.path(metadata_dir, roi_info_name)
roi_df <- read_csv(roi_info_path)
roi_df <- roi_df[roi_df$ROI_Location != "AdjacentNormal",]
# roi_df$ROI_Location[roi_df$ROI_Location=="DistantNormal"] <- "Normal"
roi_df <- roi_df[roi_df$ROI_ID %in% cell_unique_rois, ]
roi_lst <- roi_df$ROI_ID

roi_stages <- character(0)
roi_cluster_vals <- numeric(0)
for (ind in 1:length(roi_lst)) {
  roi_name <- roi_lst[ind]
  cell_roi_indices <- which(cell_rois == roi_name)
  cell_roi_locs <- cell_locs[cell_roi_indices, ]
  if (length(cell_roi_locs) < 5)
    next
  # collect information
  roi_stages <- c(roi_stages, roi_df$ROI_Diag[ind])
  roi_cluster_vals <- c(roi_cluster_vals, hopkins(cell_roi_locs))
}


# stage_cluster_exprs <- c()
# unique_stages <- c("Normal", "AAH", "AIS", "MIA", "ADC")
# for (cur_stage in unique_stages) {
#   roi_stage_indices <- which(roi_stages == cur_stage)
#   stage_cluster_exprs[cur_stage] <- roi_cluster_vals[roi_stage_indices]
# }


roi_cluster_df <- data.frame(cluster_vals=roi_cluster_vals, Stage=roi_stages)
roi_cluster_df$Stage <- factor(roi_cluster_df$Stage, levels = c("Normal", "AAH", "AIS", "MIA", "ADC"))
ggplot(roi_cluster_df, aes(x = Stage, y = cluster_vals, colour = Stage)) + 
  geom_boxplot() +
  # geom_violin() + 
  geom_jitter(width = 0.1) + 
  ylim(0.0, 1.0) +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))