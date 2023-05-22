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
cell_type_interaction_name <- paste0("DelaunayInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_delaunay_dir, cell_type_interaction_name)
cell_spatial_path <- file.path(cell_type_interaction_path)
load(cell_spatial_path)

ct_palette <- c(
    "#CAB2D6", # B-Cell
    "#6A3D9A", # CD4-T-Cell
    "#1F77B4", # CD8-T-Cell
    "#A5CEE2", # Dendritic-Cell
    "#FF7F00", # Endothelial-Cell
    "#E21A1C", # Epithelial-Cell
    "#FDBE6F", # Fibroblast
    "#33A02B", # Macrophage
    "#B2DF8A", # MDSC
    "#FEFF99", # Monocyte
    "#663300", # Neutrophil
    "#E658A0", # NK-Cell
    "#666600", # Proliferating-Cell
    "#333300", # T-Reg-Cell
    "#D6D6D6"  # Undefined
)

img_file_path <- file.path(data_root_dir, "LungROIProcessing", "Steinbock", "images.csv")
img_data_df <- read_csv(img_file_path)
num_img <- nrow(img_data_df)

ct_delaunay_dir <- file.path(data_root_dir, "Results", "VisCT-Delaunay")
if (!file.exists(ct_delaunay_dir))
    dir.create(ct_delaunay_dir)

# # test roi
# # Normal - 2323-5F-ROI002
# # AAH - 2513-1G-ROI005
# # AIS - H18-0255-6-ROI015
# # MIA - 2513-1C-ROI012
# # ADC - H18-0331-11-ROI003
# test_roi_name <- "2323-5F-ROI002"
# # cell connection visualization
# plotSpatial(spe[, spe$sample_id == test_roi_name],
#             img_id = "sample_id",
#             node_color_by = "celltype",
#             node_size_fix = 1.0,
#             draw_edges = TRUE,
#             colPairName = "delaunay_interaction_graph",
#             edge_width_fix = 0.1,
#             nodes_first = FALSE,
#             directed = FALSE,
#             edge_color_fix = "black") + 
#     scale_color_manual(values = ct_palette) 


for (ind in 1:num_img) {
    img_fullname <- img_data_df$image[ind]
    img_name <- substr(img_fullname, 1, nchar(img_fullname) - 5)
    ct_delaunay_path <- file.path(ct_delaunay_dir, paste0(img_name, ".pdf"))
    # # test roi
    # test_roi_name <- "2166-1B-ROI009"
    # cell connection visualization
    plotSpatial(spe[, spe$sample_id == img_name],
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
