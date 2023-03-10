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
library(knitr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

set.seed(1234)
spe <- aggregateNeighbors(spe, colPairName = "neighborhood", aggregate_by = "metadata", count_by = "celltype")
cn_kmeans <- kmeans(spe$aggregatedNeighbors, centers = 10)
spe$cn_celltypes <- as.factor(cn_kmeans$cluster)

# # save cn
# cell_df <- data.frame(cell_id = rownames(spe@colData), 
#                       cell_type = spe$celltype,
#                       cell_cn = spe$cn_celltypes)
# cell_type_path <- file.path(phenotype_dir, "cell_cn_types.csv")
# write.csv(cell_df, cell_type_path, row.names=FALSE)


# # Cellular neighborhood analysis heatmap
# for_plot <- prop.table(table(spe$cn_celltypes, spe$celltype), margin = 1)
# pheatmap(for_plot, color = colorRampPalette(c("dark blue", "white", "dark red"))(10), scale = "column")

# # save the spe with cellular neighborhood information
# spe_cn_path <- file.path(data_root_dir, "SpatialAnalysis", "spe_cn10.rds")
# saveRDS(spe, spe_cn_path)


# plotting image cell type information
numCN <- length(unique(spe$cn_celltypes))
colCN <- setNames(colorRampPalette(brewer.pal(10, "Set3"))(numCN), sort(unique(spe$cn_celltypes)))
colCN["2"] <- "#E41A1C"
colCN["7"] <- "#33A02C"
colCN["9"] <- "#A65628"
colCN["10"] <- "#6A3D9A"


img_file_path <- file.path(data_root_dir, "LungROIProcessing", "Steinbock", "images.csv")
img_data_df <- read_csv(img_file_path)
num_img <- nrow(img_data_df)

cn_vis_dir <- file.path(data_root_dir, "Results", "CellularNeighbor10")
if (!file.exists(cn_vis_dir))
    dir.create(cn_vis_dir)


for (ind in 1:num_img) {
    img_fullname <- img_data_df$image[ind]
    img_name <- substr(img_fullname, 1, nchar(img_fullname) - 5)
    vis_cn_path <- file.path(cn_vis_dir, paste0(img_name, ".pdf"))
    # plot
    plotSpatial(spe[, spe$sample_id == img_name],
                node_color_by = "cn_celltypes", img_id = "sample_id", node_size_fix = 1.5) + 
        scale_color_manual(values = colCN, limits = force) + 
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title=element_blank(),
              legend.position="none") 
    #     scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) 
    # save
    vis_w <- img_data_df$height_px[ind] / 100
    vis_h <- img_data_df$width_px[ind] / 100
    ggsave(filename = vis_cn_path, device='pdf', width=vis_w, height=vis_h, dpi=300)
    while (!is.null(dev.list()))
        dev.off()
}
