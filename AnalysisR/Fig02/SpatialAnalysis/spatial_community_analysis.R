library(SpatialExperiment)
library(imcRtools)
library(BiocParallel)
library(cytomapper)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(viridis)
library(pheatmap)
library(scales)
library(igraph)
library(ggraph)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)


# Spatial community detection - tumor
tumor_spe <- spe[,spe$celltype == "Epithelial-Cell"]
gr <- graph_from_data_frame(as.data.frame(colPair(tumor_spe, "neighborhood")), directed = FALSE, 
                            vertices = data.frame(index = seq_len(ncol(tumor_spe))))
cl_comm <- cluster_louvain(gr)
comm_tumor <- paste0("Tumor_", membership(cl_comm))
comm_tumor[membership(cl_comm) %in% which(sizes(cl_comm) < 10)] <- NA
names(comm_tumor) <- colnames(tumor_spe)

# Spatial community detection - non-tumor
stroma_spe <- spe[,spe$celltype != "Epithelial-Cell"]
gr <- graph_from_data_frame(as.data.frame(colPair(stroma_spe, "neighborhood")), directed = FALSE, 
                            vertices = data.frame(index = seq_len(ncol(stroma_spe))))
cl_comm <- cluster_louvain(gr)
comm_stroma <- paste0("Stroma_", membership(cl_comm))
comm_stroma[membership(cl_comm) %in% which(sizes(cl_comm) < 10)] <- NA
names(comm_stroma) <- colnames(stroma_spe)
comm <- c(comm_tumor, comm_stroma)
spe$spatial_community <- comm[colnames(spe)]


continuity_dir <- file.path(data_root_dir, "NatureFigures", "Fig02", "ROIsContinuity")
continuity_file_path <-file.path(continuity_dir, "ROIsContinuty.xlsx")
continuity_data <- read_excel(continuity_file_path)
continuity_roi_lst <- continuity_data$ROI_ID

for (continuity_roi in continuity_roi_lst) {
    file_continuity_path <- file.path(continuity_dir, "EpithelialCommunity", paste0(continuity_roi, ".pdf"))
    roi_spe <- spe[, spe$sample_id == continuity_roi]
    
    plotSpatial(roi_spe[, roi_spe$celltype == "Epithelial-Cell"], 
                node_color_by = "spatial_community", 
                img_id = "sample_id", node_size_fix = 1.0) +
                theme(legend.position = "none")
    ggsave(filename = file_continuity_path, device='pdf', width=10, height=9, dpi=300)
    while (!is.null(dev.list()))  
        dev.off()    
}


