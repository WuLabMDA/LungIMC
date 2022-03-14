# batch effect correction visualization
vis_correction_dir <- file.path(data_root, "BatchCorrection", "VisCorrection")
if (!dir.exists(vis_correction_dir))
    dir.create(vis_correction_dir, recursive = TRUE)

# Create UMAPs
sam <- sample(1:nrow(uncorrected), 50000)
sample_uncorrected <- uncorrected[sam, ]
sample_corrected <- corrected[sam, ]
plot1 <- plot_dimred(sample_uncorrected, "Uncorrected", type = "umap", plot = "batch", markers = markers)
plot2 <- plot_dimred(sample_corrected, "Corrected", type = "umap", plot = "batch", markers = markers)
# Umap plots
plot_save_two(plot1, plot2, file.path(vis_correction_dir, "umap.png"))
# Density plots
plot_density(sample_uncorrected, sample_corrected, markers = markers,
             filename = file.path(vis_correction_dir, "density.png"), 
             y = "batch", ncol = 6)