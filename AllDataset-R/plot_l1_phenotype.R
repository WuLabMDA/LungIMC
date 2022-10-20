library(SpatialExperiment)
library(patchwork)
library(Seurat)
library(SeuratObject)
library(dittoSeq)
library(RColorBrewer)
library(viridis)

set.seed(1234)
data_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockAll"
phenotype_path <- file.path(data_root_dir, "som_l1_spe.rds")
spe <- readRDS(phenotype_path)

fig_dir <- file.path(data_root_dir, "L1Figs")
if (!dir.exists(fig_dir)){
    dir.create(fig_dir)
}

#  UMAP embedding
umap_plot_path <- file.path(fig_dir, "som_l1_umap.png")
png(file=umap_plot_path, width=1200, height=1000, units = "px")
cur_cells <- sample(seq_len(ncol(spe)), 100000)
dittoDimPlot(spe[, cur_cells], var = "l1_clusters", reduction.use = "UMAP", size = 0.1, do.label = TRUE) +
    ggtitle("SOM clusters expression on UMAP")
dev.off()

# dot plot
seurat_obj <- as.Seurat(spe, counts = "counts", data = "exprs")
seurat_obj <- AddMetaData(seurat_obj, as.data.frame(colData(spe)))
dot_plot_path <- file.path(fig_dir, "som_l1_dotplot.png")
png(file=dot_plot_path, width=1200, height=1000, units = "px")
l1_feas <- rownames(spe)[rowData(spe)$l1_markers]
DotPlot(seurat_obj, features = l1_feas, group.by = 'l1_clusters', cols = c("#deebf7", "#08519c"), dot.scale = 9) + 
    RotatedAxis() + theme(legend.position = "bottom") + labs(x = NULL, y = NULL)
dev.off()
