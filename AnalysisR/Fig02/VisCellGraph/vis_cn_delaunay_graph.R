library(imcRtools)
library(SpatialExperiment)
library(stringr)
library(RColorBrewer)
library(dittoSeq)
library(ggplot2)
library(viridis)
library(randomcoloR)
library(tidyverse)
library(readxl)
library(knitr)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

celltype_delaunay_dir <- file.path(phenotype_dir, "DelaunayInteraction")
threshold_val <- 50
num_cn <- 8
spe_cn_name <- paste0("Delaunay", threshold_val, "-CN", num_cn, ".rds")
spe_cn_path <- file.path(celltype_delaunay_dir, spe_cn_name)
spe <- readRDS(spe_cn_path)

cn_palette <- c(
    "#D6D6D6", # Undefined-CN1
    "#E21A1C", # Epithelial1-CN2
    "#666600", # Proliferating-CN3
    "#FF7F00", # Epithelial2-CN4
    "#FF20A9", # Endothelial-CN5
    "#FDBE6F", # Fibroblast-CN6
    "#33A02B", # Macrophage-CN7
    "#B2DF8A" # PanImmune-CN8
)

img_file_path <- file.path(data_root_dir, "LungROIProcessing", "Steinbock", "images.csv")
img_data_df <- read_csv(img_file_path)
num_img <- nrow(img_data_df)

cn_delaunay_dir <- file.path(data_root_dir, "Results", "VisCN-Delaunay")
if (!file.exists(cn_delaunay_dir))
    dir.create(cn_delaunay_dir)

for (ind in 1:num_img) {
    img_fullname <- img_data_df$image[ind]
    img_name <- substr(img_fullname, 1, nchar(img_fullname) - 5)
    cn_delaunay_path <- file.path(cn_delaunay_dir, paste0(img_name, ".pdf"))
    # cell connection visualization
    plotSpatial(spe[, spe$sample_id == img_name],
                img_id = "sample_id",
                node_color_by = "cn_celltypes",
                node_size_fix = 1.0,
                draw_edges = TRUE,
                colPairName = "delaunay_interaction_graph",
                edge_width_fix = 0.1,
                nodes_first = FALSE,
                directed = FALSE,
                edge_color_fix = "black") +
        scale_color_manual(values = cn_palette) +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title=element_blank(),
              legend.position="none")
    # save
    vis_w <- img_data_df$height_px[ind] / 100
    vis_h <- img_data_df$width_px[ind] / 100
    ggsave(filename = cn_delaunay_path, device='pdf', width=vis_w, height=vis_h, dpi=300)
    while (!is.null(dev.list()))
        dev.off()
}