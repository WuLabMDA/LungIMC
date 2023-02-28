library(SpatialExperiment)
library(imcRtools)
library(cytomapper)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(viridis)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(readxl)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

set.seed(1234)
spe <- aggregateNeighbors(spe, colPairName = "neighborhood", 
                          aggregate_by = "metadata", count_by = "celltype")
cn_kmeans <- kmeans(spe$aggregatedNeighbors, centers = 24)
spe$cn_celltypes <- as.factor(cn_kmeans$cluster)

# # save the spe with cellular neighborhood information
# spe_cn_path <- file.path(data_root_dir, "SpatialAnalysis", "spe_cn24.rds")
# saveRDS(spe, spe_cn_path)

# Cellular neighborhood analysis heatmap
for_plot <- prop.table(table(spe$cn_celltypes, spe$celltype), margin = 1)
pheatmap(for_plot, color = colorRampPalette(c("dark blue", "white", "dark red"))(100), scale = "column")


# plotting image cell type information
numCN <- length(unique(spe$cn_celltypes))
colCN <- setNames(colorRampPalette(brewer.pal(9, "Paired"))(numCN), sort(unique(spe$cn_celltypes)))
continuity_dir <- file.path(data_root_dir, "NatureFigures", "Fig02", "ROIsContinuity")
continuity_file_path <-file.path(continuity_dir, "ROIsContinuty.xlsx")
continuity_data <- read_excel(continuity_file_path)

continuity_roi_lst <- continuity_data$ROI_ID
for (continuity_roi in continuity_roi_lst) {
    file_continuity_path <- file.path(continuity_dir, paste0(continuity_roi, ".pdf"))
    plotSpatial(spe[, spe$sample_id == continuity_roi], 
            node_color_by = "cn_celltypes", 
            img_id = "sample_id", 
            node_size_fix = 1.0) +
    scale_color_manual(values = colCN, limits = force)
    ggsave(filename = file_continuity_path, device='pdf', width=10, height=8, dpi=300)
    while (!is.null(dev.list()))  
        dev.off()
}
