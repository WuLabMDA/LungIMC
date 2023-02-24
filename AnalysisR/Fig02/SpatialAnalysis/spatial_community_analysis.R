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


## extract the front 6
unique_roi_lst <-unique(spe$sample_id)
uniuqe_roi_num <- length(unique_roi_lst)
cell_id_lst <- rownames(colData(spe))
roi_num <- 6
ttl_cell_num <- 0
for (ir in 1:roi_num) {
    cur_roi = unique_roi_lst[ir]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_cell_num <- length(cell_indices)
    ttl_cell_num <- ttl_cell_num + roi_cell_num
}
display_spe <- spe[, 1:ttl_cell_num]

# display
plotSpatial(display_spe[, display_spe$celltype == "Epithelial-Cell"], 
            node_color_by = "spatial_community", 
            img_id = "sample_id", node_size_fix = 1.0)



