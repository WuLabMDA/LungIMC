# Batch correct
corrected <- batch_correct(df = uncorrected, covar = "condition", markers = imc_markers,
                           norm_method = "scale", rlen = 10, seed = 1234)
# Re-run clustering on corrected data
labels <- corrected %>% create_som(markers = imc_markers, rlen = 10)
uncorrected$label <- corrected$label <- labels

emd <- evaluate_emd(uncorrected, corrected, cell_col = "label")
mad <- evaluate_mad(uncorrected, corrected, cell_col = "label")