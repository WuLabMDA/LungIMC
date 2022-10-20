library(SpatialExperiment)
library(dittoSeq)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(Seurat)
library(SeuratObject)

set.seed(1234)
data_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockAll"

phenotype_path <- file.path(data_root_dir, "som_raw_spe.rds")
spe <- readRDS(phenotype_path)
seurat_obj <- as.Seurat(spe, counts = "counts", data = "exprs")
seurat_obj <- AddMetaData(seurat_obj, as.data.frame(colData(spe)))

fig_dir <- file.path(data_root_dir, "Figs")
if (!dir.exists(fig_dir)){
    dir.create(fig_dir)
}

# Draw dot plot on Raw 
dotplot_name <- "som_raw_dotplot"
plot_path <- file.path(fig_dir, paste(dotplot_name, ".png", sep=""))
png(file=plot_path, width=1200, height=1000, units = "px")
interested_feas <- rownames(spe)[rowData(spe)$use_channel]
DotPlot(seurat_obj, features = interested_feas, group.by = 'som_clusters', cols = c("#deebf7", "#08519c"), dot.scale = 9) + 
    RotatedAxis() + theme(legend.position = "bottom") + labs(x = NULL, y = NULL)
dev.off()
