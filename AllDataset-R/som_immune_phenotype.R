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
library(dittoSeq)
library(viridis)
library(patchwork)
library(Seurat)
library(SeuratObject)
set.seed(1234)

# load all data clustering
data_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockAll"
all_phenotype_path <- file.path(data_root_dir, "som_lineage_spe.rds")
all_spe <- readRDS(all_phenotype_path)

# obtain cluster information
all_cluster_ids <- all_spe$som_clusters 
# filtering immune clusters [3, 6, 8, 10, 11]
immnue_clusters <- c(3, 6, 8, 10, 11)
cell_indices <- which(all_cluster_ids %in% immnue_clusters)
immnue_spe <- all_spe[, cell_indices]

# Define lineage and function marker list
# immnue_strings <- "CD45|CD3e|CD4|FoxP3|CD8a|CD19|CD94|CD11b|CD11c|CD14|MPO|CD68|CD33|HLADR|CD45RO|Ki67|CD163|ICOS|GranzymeB|TIM3|LAG3|VISTA|PD1|PDL1|TIGIT|IDO1|B7H3|CTLA4|CD73"
immnue_strings <- "CD45|CD3e|CD4|FoxP3|CD8a|CD19|CD94|CD11b|CD11c|CD14|MPO|CD68|CD33|HLADR|CD45RO"
rowData(immnue_spe)$immnue_markers <- grepl(immnue_strings, rownames(immnue_spe))

# Run FlowSOM and ConsensusClusterPlus clustering
immnue_spe <- cluster(immnue_spe, features = rownames(immnue_spe)[rowData(immnue_spe)$immnue_markers], maxK = 50, seed = 1234)
# # Assess cluster stability
# delta_area(immnue_spe)

immnue_spe$immnue_clusters <- cluster_ids(immnue_spe, "meta12")
# # save som clustered results
# immnue_phenotype_path <- file.path(data_root_dir, "som_immnue_spe.rds")
# saveRDS(immnue_spe, immnue_phenotype_path)
# create figure folders
fig_dir <- file.path(data_root_dir, "ImmnueFigs")
if (!dir.exists(fig_dir)){
    dir.create(fig_dir)
}
# Convert to Seurat
seurat_obj <- as.Seurat(immnue_spe, counts = "counts", data = "exprs")
seurat_obj <- AddMetaData(seurat_obj, as.data.frame(colData(immnue_spe)))
# Draw dot plot on Raw 
dotplot_name <- "som_lineage_dotplot"
plot_path <- file.path(fig_dir, paste(dotplot_name, ".png", sep=""))
png(file=plot_path, width=1200, height=1000, units = "px")
interested_feas <- rownames(immnue_spe)[rowData(immnue_spe)$immnue_markers]
DotPlot(seurat_obj, features = interested_feas, group.by = 'immnue_clusters', cols = c("#deebf7", "#08519c"), dot.scale = 9) + 
    RotatedAxis() + theme(legend.position = "bottom") + labs(x = NULL, y = NULL)
dev.off()
#  UMAP embedding
cluster_plot_name <- "som_lineage_umap"
plot_path <- file.path(fig_dir, paste(cluster_plot_name, ".png", sep=""))
png(file=plot_path, width=1200, height=1000, units = "px")
cur_cells <- sample(seq_len(ncol(immnue_spe)), 100000)
dittoDimPlot(immnue_spe[, cur_cells], var = "immnue_clusters",  reduction.use = "UMAP", size = 0.2, do.label = TRUE) +
    ggtitle("SOM clusters expression on UMAP")
dev.off()



