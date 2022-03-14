# batch effect correction visualization
vis_correction_dir <- file.path(data_root, "BatchCorrection", "NoControl")
if (!dir.exists(vis_correction_dir))
    dir.create(vis_correction_dir, recursive = TRUE)

# save objects
no_control_meta_path <- file.path(vis_correction_dir, "meta.RData")
save(fea_filenames, markers, roi_nrows, file = no_control_meta_path)
raw_fea_path <- file.path(vis_correction_dir, "RawFeas.RData")
save(roi_feas, file = raw_fea_path)
transform_fea_path <- file.path(vis_correction_dir, "TransformedFeas.RData")
save(uncorrected, file = transform_fea_path)
correct_fea_path <- file.path(vis_correction_dir, "CorrectedFeas.RData")
save(corrected, file = correct_fea_path)