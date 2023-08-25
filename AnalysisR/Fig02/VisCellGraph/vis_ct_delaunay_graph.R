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
library(comprehenr)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

celltype_delaunay_dir <- file.path(phenotype_dir, "DelaunayInteraction")
threshold_val <- 50
cell_type_interaction_name <- paste0("DelaunayInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_delaunay_dir, cell_type_interaction_name)
cell_spatial_path <- file.path(cell_type_interaction_path)
load(cell_spatial_path)


ct_palette <- c(
    "#9467BD", # B-Cell
    "#2CA02C", # CD4-T-Cell
    "#98DF8A", # CD8-T-Cell
    "#FDBE6F", # Dendritic-Cell
    "#FF9F00", # Endothelial-Cell
    "#6CCEE2", # Epithelial-Cell
    "#E387B4", # Fibroblast
    "#8C754B", # Macrophage
    "#B2DF8A", # MDSC
    "#C49C94", # Monocyte
    "#7F7F7F", # Neutrophil
    "#DBDB8D", # NK-Cell
    "#769AE1", # Proliferating-Cell
    "#FF9896", # T-Reg-Cell
    "#D6D6D6"  # Undefined
)


# map_ct_palette <- c("B-Cell"="#9467BD", "CD4-T-Cell"="#2CA02C", "CD8-T-Cell"="#98DF8A",
#                     "Dendritic-Cell"="#FDBE6F", "Endothelial-Cell"="#FF9F00", "Epithelial-Cell"="#6CCEE2",
#                     "Fibroblast"="#E387B4", "Macrophage"="#8C754B", "MDSC"="#B2DF8A",
#                     "Monocyte"="#C49C94", "Neutrophil"="#7F7F7F", "NK-Cell"="#DBDB8D",
#                     "Proliferating-Cell"="#769AE1", "T-Reg-Cell"="#FF9896", "Undefined"="#D6D6D6")

img_file_path <- file.path(data_root_dir, "LungROIProcessing", "Steinbock", "images.csv")
img_data_df <- read_csv(img_file_path)
num_img <- nrow(img_data_df)

ct_delaunay_dir <- file.path(data_root_dir, "Results", "VisCT-Delaunay")
if (!file.exists(ct_delaunay_dir))
    dir.create(ct_delaunay_dir)


for (ind in 1:num_img) {
    img_fullname <- img_data_df$image[ind]
    img_name <- substr(img_fullname, 1, nchar(img_fullname) - 5)
    ct_delaunay_path <- file.path(ct_delaunay_dir, paste0(img_name, ".pdf"))
    # # test roi
    # test_roi_name <- "2166-1B-ROI009"
    # cell connection visualization
    roi_object <- spe[, spe$sample_id == img_name]
    # roi_celltypes <- unique(object$celltype)
    # ct_palette <- to_vec(for (ct in roi_celltypes) map_ct_palette[ct])
    plotSpatial(roi_object,
                img_id = "sample_id",
                node_color_by = "celltype",
                node_size_fix = 1.0,
                draw_edges = TRUE,
                colPairName = "delaunay_interaction_graph",
                edge_width_fix = 0.1,
                nodes_first = FALSE,
                directed = FALSE,
                edge_color_fix = "black") +
        scale_color_manual(values = ct_palette) +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title=element_blank(),
              legend.position="none")
    # save
    vis_w <- img_data_df$height_px[ind] / 100
    vis_h <- img_data_df$width_px[ind] / 100
    ggsave(filename = ct_delaunay_path, device='pdf', width=vis_w, height=vis_h, dpi=300)
    while (!is.null(dev.list()))
        dev.off()
}
