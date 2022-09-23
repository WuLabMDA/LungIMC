library(SpatialExperiment)
library(dittoSeq)
library(RColorBrewer)
library(viridis)

set.seed(1234)
normal_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockAll"

fig_dir <- file.path(normal_root_dir, "Figs")
if (!dir.exists(fig_dir)){
    dir.create(fig_dir)
}

phenotype_path <- file.path(normal_root_dir, "som_raw_spe.rds")
spe <- readRDS(phenotype_path)

#  UMAP embedding
cluster_plot_name <- "som_raw_umap"
plot_path <- file.path(fig_dir, paste(cluster_plot_name, ".png", sep=""))
png(file=plot_path, width=1200, height=1000, units = "px")
dittoDimPlot(spe, var = "som_clusters",  reduction.use = "UMAP", size = 0.2, do.label = TRUE) +
    ggtitle("SOM clusters expression on UMAP")
dev.off()
