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

# load the spe
threshold_val <- 50
num_cn <- 10
spe_cn_name <- paste0("Delaunay", threshold_val, "-CN", num_cn, ".rds")
spe_cn_path <- file.path(data_root_dir, "CellPhenotyping", spe_cn_name)
spe <- readRDS(spe_cn_path)

# plotting image cell type information
colCN <- setNames(colorRampPalette(brewer.pal(10, "Set3"))(num_cn), sort(unique(spe$cn_celltypes)))
colCN["2"] <- "#E41A1C"
colCN["7"] <- "#33A02C"
colCN["9"] <- "#A65628"
colCN["10"] <- "#6A3D9A"


img_file_path <- file.path(data_root_dir, "LungROIProcessing", "Steinbock", "images.csv")
img_data_df <- read_csv(img_file_path)
num_img <- nrow(img_data_df)

cn_vis_dir <- file.path(data_root_dir, "Results", "Vis-CN10-Delaunay")
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
    # save
    vis_w <- img_data_df$height_px[ind] / 100
    vis_h <- img_data_df$width_px[ind] / 100
    ggsave(filename = vis_cn_path, device='pdf', width=vis_w, height=vis_h, dpi=300)
    while (!is.null(dev.list()))
        dev.off()
}
