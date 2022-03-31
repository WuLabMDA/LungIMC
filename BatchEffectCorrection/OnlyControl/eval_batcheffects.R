# set batch effects evaluation directory
batch_effects_dir <- file.path(data_root, "EvalBatchEffects", "Control")
if (!dir.exists(batch_effects_dir))
    dir.create(batch_effects_dir, recursive = TRUE)

# checking for batch effects
detect_batch_effect(uncorrected, out_dir = batch_effects_dir, downsample = NULL,
                    norm_method = "scale", xdim = 8, ydim = 8,
                    markers = markers, batch_col = "batch", seed = 1234)

# quick detection of batch effects.
detect_batch_effect_express(uncorrected, out_dir = batch_effects_dir,
                            batch_col = "batch", downsample = NULL, seed = 1234)
