# Batch correct
corrected <- batch_correct(df = uncorrected, covar = "condition", markers = markers,
                           norm_method = "scale", rlen = 10, seed = 1234)
# Re-run clustering on corrected data
labels <- corrected %>% create_som(markers = markers, rlen = 10)
uncorrected$label <- corrected$label <- labels

## save corrected features
rdata_dir <- file.path(data_root, "RData")
if (!dir.exists(rdata_dir))
    dir.create(rdata_dir, recursive = TRUE)
CorrectedFeas <- corrected %>% data.frame()
correct_fea_path <- file.path(rdata_dir, "CorrectedFeas.csv")
write.csv(CorrectedFeas, correct_fea_path, row.names = FALSE)


# batch effect correction visualization
vis_correction_dir <- file.path(data_root, "NoControl")
if (!dir.exists(vis_correction_dir))
    dir.create(vis_correction_dir, recursive = TRUE)

# Evaluate MAD
mad <- evaluate_mad(uncorrected, corrected, cell_col = "label")
message(paste("MAD score is:", mad$score))
# Evaluate EMD
emd <- evaluate_emd(uncorrected, corrected, cell_col = "label")
message(paste("EMD Reduction is:", emd$reduction))
png(file=file.path(vis_correction_dir, "emd_reduction.png"), width=1600, height=1000, res=180)
cowplot::plot_grid(emd$violin, emd$scatterplot)
dev.off()

# Create UMAPs
sam <- sample(1:nrow(uncorrected), 100000)
sample_uncorrected <- uncorrected[sam, ]
sample_corrected <- corrected[sam, ]
plot1 <- plot_dimred(sample_uncorrected, "Uncorrected", type = "umap", plot = "batch", markers = markers)
plot2 <- plot_dimred(sample_corrected, "Corrected", type = "umap", plot = "batch", markers = markers)
# Umap plots
plot_save_two(plot1, plot2, file.path(vis_correction_dir, "umap.png"))
# Density plots
plot_density(sample_uncorrected, sample_corrected, markers = markers, y = "batch", ncol = 6,
             filename = file.path(vis_correction_dir, "density.png"))
