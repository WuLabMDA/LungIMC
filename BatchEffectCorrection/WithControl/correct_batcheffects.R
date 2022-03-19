# Batch correct
corrected <- batch_correct(df = uncorrected, covar = "condition",
                           anchor = "anchor", markers = markers,
                           norm_method = "scale", rlen = 10, seed = 1234)
uncorrected <- uncorrected %>% filter(anchor != "Replicates")
corrected <- corrected %>% filter(anchor != "Replicates")
# run clustering on corrected data
labels <- corrected %>% create_som(markers = markers, rlen = 10)
uncorrected$label <- corrected$label <- labels

## save corrected features
rdata_dir <- file.path(data_root, "RData")
if (!dir.exists(rdata_dir))
    dir.create(rdata_dir, recursive = TRUE)
CorrectedFeas <- corrected %>% data.frame()
correct_fea_path <- file.path(rdata_dir, "ControlCorrectedFeas.csv")
write.csv(CorrectedFeas, correct_fea_path, row.names = FALSE)