# load Panel metadata
panel_file_path <- file.path(data_root, "Metadata", "PanelsIMC.csv")
markers <- read.csv(panel_file_path) %>% pull(Marker)


# load ROI metadata
roi_meta_path <- file.path(data_root, "Metadata", "NoControlMetaIMC.csv")
roi_metadata <- suppressMessages(readr::read_csv(roi_meta_path))
fea_filenames <- roi_metadata$Filename


# load features
file_num <- length(fea_filenames)
roi_fea_dir <- file.path(data_root, "NoControl", "CellFeas")
np <- import("numpy")
p_len <- 2
p_num <- 35
roi_fea_list <- vector(mode = "list", length = file_num)
roi_nrows <- vector(length = file_num)
for (ind in 1:file_num) {
    cur_fea_path = file.path(roi_fea_dir, paste0(fea_filenames[ind], ".npy"))
    cur_fea_mat <- tibble::as_tibble(np$load(cur_fea_path)[, (p_len+1):(p_len+p_num)])
    roi_fea_list[[ind]] <- cur_fea_mat
    roi_nrows[ind] <- nrow(cur_fea_mat)
}

# convert to feature matrix and add markers
roi_feas <- bind_rows(roi_fea_list)
colnames(roi_feas) <- markers
# add cell ids
roi_feas$id <- 1:sum(roi_nrows)
# add ROI ids
roi_feas$sample <- fea_filenames %>% rep(roi_nrows)
# add batch info
batch_ids <- roi_metadata[["Batch"]][match(fea_filenames, roi_metadata[["Filename"]])] %>% as.factor() %>% rep(roi_nrows)
roi_feas$batch <- batch_ids
# add condition info
condition_ids <- roi_metadata[["Condition"]][match(fea_filenames, roi_metadata[["Filename"]])] %>% as.factor() %>% rep(roi_nrows)
roi_feas$condition <- condition_ids

# transform data
uncorrected <- transform_asinh(df = roi_feas, markers = markers, cofactor = 5, .keep = TRUE)

# save RData
rdata_dir <- file.path(data_root, "RData")
if (!dir.exists(rdata_dir))
    dir.create(rdata_dir, recursive = TRUE)
## save meta data
no_control_meta_path <- file.path(rdata_dir, "StudyMeta.RData")
save(markers, fea_filenames, roi_nrows, file = no_control_meta_path)
# ## save raw features
# raw_labels <- roi_feas %>% create_som(markers = markers, rlen = 10)
# roi_feas$label <- raw_labels
# RawFeas <- roi_feas %>% data.frame()
# raw_fea_path <- file.path(rdata_dir, "RawFeas.csv")
# write.csv(RawFeas, raw_fea_path, row.names = FALSE)
## save transformed features
transform_labels <- uncorrected %>% create_som(markers = markers, rlen = 10)
uncorrected$label <- transform_labels
TransformFeas <- uncorrected %>% data.frame()
transform_fea_path <- file.path(rdata_dir, "TransformFeas.csv")
write.csv(TransformFeas, transform_fea_path, row.names = FALSE)
