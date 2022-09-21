library(cowplot)
library(dittoSeq)
library(viridis)
library(scater)

normal_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockDistantNormal"
seurat_spe_path = file.path(normal_root_dir, "seurat_spe.rds")
spe <- readRDS("./seurat_spe.rds")

fig_dir <- file.path(normal_root_dir, "Figs")
if (!dir.exists(fig_dir)){
    dir.create(fig_dir)
}

# Visualize correction
plot_name <- "seurat_correction"
plot_path <- file.path(fig_dir, paste(plot_name, ".pdf", sep=""))
pdf(file=plot_path)
p1 <- dittoDimPlot(spe, var = "stain_id", reduction.use = "UMAP", size = 0.2) +
    scale_color_manual(values = metadata(spe)$color_vectors$Stain) +
    ggtitle("UMAP Before Correction")
p2 <- dittoDimPlot(spe, var = "stain_id", reduction.use = "UMAP_seurat", size = 0.2) +
    scale_color_manual(values = metadata(spe)$color_vectors$Stain) +
    ggtitle("UMAP After Correction")
plot_grid(p1, p2)
dev.off()