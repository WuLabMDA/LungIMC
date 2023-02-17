library(imcRtools)
library(SpatialExperiment)
library(stringr)
library(RColorBrewer)
library(dittoSeq)
library(ggplot2)
library(viridis)
library(pheatmap)
library(igraph)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

## load SpatialExperiment data
steinbock_spe_name <- "steinbock_spe_2023_02_10"
spe_file_path <- file.path(phenotype_dir, paste0(steinbock_spe_name, ".rds"))
load(spe_file_path)
data$celltype <- data@metadata$celltypes # reorganize cell type information
spe <- data


# Spatial community detection - tumor
tumor_spe <- spe[,spe$celltype == "Epithelial"]
gr <- graph_from_data_frame(as.data.frame(colPair(tumor_spe, "neighborhood")), 
                            directed = FALSE, 
                            vertices = data.frame(index = seq_len(ncol(tumor_spe))))
cl_comm <- cluster_louvain(gr)
comm_tumor <- paste0("Epithelial_", membership(cl_comm))
comm_tumor[membership(cl_comm) %in% which(sizes(cl_comm) < 10)] <- NA
names(comm_tumor) <- colnames(tumor_spe)

# Spatial community detection - non-tumor
stroma_spe <- spe[,spe$celltype != "Epithelial"]
gr <- graph_from_data_frame(as.data.frame(colPair(stroma_spe, "neighborhood")), 
                            directed = FALSE, 
                            vertices = data.frame(index = seq_len(ncol(stroma_spe))))
cl_comm <- cluster_louvain(gr)
comm_stroma <- paste0("Stroma_", membership(cl_comm))
comm_stroma[membership(cl_comm) %in% which(sizes(cl_comm) < 10)] <- NA
names(comm_stroma) <- colnames(stroma_spe)
comm <- c(comm_tumor, comm_stroma)
spe$spatial_community <- comm[colnames(spe)]


## display the front 10
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


plotSpatial(display_spe[, display_spe$celltype == "Epithelial"], 
            node_color_by = "spatial_community", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
    theme(legend.position = "none") +
    ggtitle("Spatial tumor communities") +
    scale_color_manual(values = rev(colors()))

# By celltypes
spe <- aggregateNeighbors(spe, colPairName = "neighborhood", 
                          aggregate_by = "metadata", count_by = "celltype")
cn_1 <- kmeans(spe$aggregatedNeighbors, centers = 10)
spe$cn_celltypes <- as.factor(cn_1$cluster)


plotSpatial(spe[, 1:ttl_cell_num], 
            node_color_by = "cn_celltypes", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
    scale_color_brewer(palette = "Set3")

for_plot <- prop.table(table(spe$cn_celltypes, spe$celltype), margin = 1)
pheatmap(for_plot, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column")
