library(imcRtools)
library(SpatialExperiment)
library(CATALYST)
library(cytomapper)
library(scran)
library(scater)
library(bluster)
library(BiocParallel)
library(kohonen)
library(ConsensusClusterPlus)
library(Rphenograph)
library(igraph)
library(dittoSeq)
library(viridis)
library(patchwork)
library(Seurat)
library(SeuratObject)
set.seed(1234)

# load all data clustering
data_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockAll"
all_phenotype_path <- file.path(data_root_dir, "som_l1_spe.rds")
all_spe <- readRDS(all_phenotype_path)

# Filtering by immune clusters
all_cluster_ids <- all_spe$l1_clusters 
l1_immnue_clusters <- c(3, 6, 8, 10, 11) # filtering immune clusters [3, 6, 8, 10, 11]
cell_indices <- which(all_cluster_ids %in% l1_immnue_clusters)
immnue_spe <- all_spe[, cell_indices]

# Define immune lineage marker list
# l2_antibodies <- "CD45|CD3e|CD4|FoxP3|CD8a|CD19|CD94|CD11b|CD11c|CD14|MPO|CD68|CD33|HLADR|CD45RO"
l2_antibodies <- "CD3e|^CD4$|FoxP3|CD8a|CD19|CD94|CD11b|CD11c|CD14|MPO|CD68|CD33"
rowData(immnue_spe)$l2_markers <- grepl(l2_antibodies, rownames(immnue_spe))

# Run FlowSOM and ConsensusClusterPlus clustering
immnue_spe <- cluster(immnue_spe, features = rownames(immnue_spe)[rowData(immnue_spe)$l2_markers], maxK = 50, seed = 1234)
# # Assess cluster stability
# delta_area(immnue_spe)
immnue_spe$l2_clusters <- cluster_ids(immnue_spe, "meta20")


# create figure folders
fig_dir <- file.path(data_root_dir, "L2Figs")
if (!dir.exists(fig_dir)){
    dir.create(fig_dir)
}

#  UMAP embedding
umap_plot_path <- file.path(fig_dir, "som_l2_umap_20.png")
png(file=umap_plot_path, width=1200, height=1000, units = "px")
cur_cells <- sample(seq_len(ncol(immnue_spe)), 100000)
dittoDimPlot(immnue_spe[, cur_cells], var = "l2_clusters",  reduction.use = "UMAP", size = 0.1, do.label = TRUE) +
    ggtitle("SOM clusters expression on UMAP")
dev.off()

# dot plot
seurat_obj <- as.Seurat(immnue_spe, counts = "counts", data = "exprs")
seurat_obj <- AddMetaData(seurat_obj, as.data.frame(colData(immnue_spe)))
dot_plot_path <- file.path(fig_dir, "som_l2_dotplot_20.png")
png(file=dot_plot_path, width=1200, height=1000, units = "px")
# l2_feas <- rownames(immnue_spe)[rowData(immnue_spe)$l2_markers]
l2_feas <- c("CD45", "CD3e", "CD4", "FoxP3", "CD8a", "CD19", "CD94", "CD11b", "CD11c", "CD14", "MPO", "CD68", "CD33", "HLADR", "CD45RO")
DotPlot(seurat_obj, features = l2_feas, group.by = 'l2_clusters', cols = c("#deebf7", "#08519c"), dot.scale = 9) + 
    RotatedAxis() + theme(legend.position = "bottom") + labs(x = NULL, y = NULL)
dev.off()


# Heatmap
heatmap_plot_path <- file.path(fig_dir, "som_l2_heatmap_20.png")
png(file=heatmap_plot_path, width=1200, height=1000, units = "px")
l2_feas <- c("CD45", "CD3e", "CD4", "FoxP3", "CD8a", "CD19", "CD94", "CD11b", "CD11c", "CD14", "MPO", "CD68", "CD33", "HLADR", "CD45RO")
dittoHeatmap(immnue_spe, genes = l2_feas, assay = "exprs",
             scaled.to.max = TRUE, heatmap.colors.max.scaled = colorRampPalette(c("white", "blue"))(50),
             heatmap.colors = viridis(100), annot.by = "l2_clusters", annot.colors = dittoColors(1)[1:length(unique(immnue_spe$l2_clusters))])
dev.off()