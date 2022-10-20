library(SpatialExperiment)
library(dittoSeq)
library(RColorBrewer)
library(viridis)

set.seed(1234)
data_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockAll"

fig_dir <- file.path(data_root_dir, "Figs")
if (!dir.exists(fig_dir)){
    dir.create(fig_dir)
}

phenotype_path <- file.path(data_root_dir, "som_lineage_spe.rds")
spe <- readRDS(phenotype_path)

#  UMAP embedding
cluster_plot_name <- "som_lineage_umap"
plot_path <- file.path(fig_dir, paste(cluster_plot_name, ".png", sep=""))
png(file=plot_path, width=1200, height=1000, units = "px")
cur_cells <- sample(seq_len(ncol(spe)), 100000)
dittoDimPlot(spe[, cur_cells], var = "som_clusters",  reduction.use = "UMAP", size = 0.2, do.label = TRUE) +
    ggtitle("SOM clusters expression on UMAP")
dev.off()
