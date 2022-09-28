library(cowplot)
library(dittoSeq)
library(viridis)

set.seed(1234)
normal_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockDistantNormal"
raw_spe_path = file.path(normal_root_dir, "steinbock_spe.rds")
spe <- readRDS(raw_spe_path)

fig_dir <- file.path(normal_root_dir, "Figs")
if (!dir.exists(fig_dir)){
    dir.create(fig_dir)
}

# Draw dot plot on Raw 
normal_smoke_name <- "raw_smoke_status"
plot_path <- file.path(fig_dir, paste(normal_smoke_name, ".png", sep=""))
png(file=plot_path, width=1200, height=1000, units = "px")
dittoDimPlot(spe, var = "smoke_status", reduction.use = "UMAP", size = 0.2) + 
    scale_color_manual(values = metadata(spe)$color_vectors$Smoke) +
    ggtitle("Patient ID on UMAP before correction")
dev.off()