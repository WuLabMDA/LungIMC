library(SpatialExperiment)
library(imcRtools)
library(BiocParallel)
library(cytomapper)
library(tidyverse)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
feature_dir <- file.path(data_root_dir, "FeatureAnalysis")

# load the spe
threshold_val <- 50
num_cn <- 10
spe_cn_name <- paste0("Delaunay", threshold_val, "-CN", num_cn, ".rds")
spe_cn_path <- file.path(phenotype_dir, spe_cn_name)
spe <- readRDS(spe_cn_path)

# extract cell information
cell_id_lst <- rownames(colData(spe))
cell_cn_lst <- spe$cn_celltypes
all_cn_lst <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# filter functional markers
function_markers = c("Ki67", "HLADR", "B2M", "CD45RO", "ICOS", "GranzymeB",
                     "TIM3", "LAG3", "VISTA", "PDL1", "PD1", "TIGIT",
                     "IDO1", "B7H3", "CTLA4", "CD73", "CD163")
cell_fun_exprs <- t(assay(spe, "exprs")[function_markers, ])

roi_lst <- unique(spe$sample_id)
function_fea_num <- length(all_cn_lst) * length(function_markers)
roi_fea_df <- data.frame(matrix(ncol = function_fea_num, nrow = length(roi_lst)))
# traverse ROI one-by-one
for (ind in 1:length(roi_lst)) {
    # extract cell index in each ROI
    cur_roi <- roi_lst[ind]
    roi_cell_indices <- startsWith(cell_id_lst, cur_roi)
    roi_cell_cns <- cell_cn_lst[roi_cell_indices]
    roi_cell_exprs <- cell_fun_exprs[roi_cell_indices, ]
    # aggreate functional marker states
    roi_state_feas <- c()
    for (cur_cn_name in all_cn_lst) {
        cur_roi_cellname_indices <- which(roi_cell_cns == cur_cn_name)
        if (length(cur_roi_cellname_indices) > 1) {
            cur_roi_cellname_exprs <- roi_cell_exprs[cur_roi_cellname_indices, ]
            cur_roi_cellname_means <- colMeans(cur_roi_cellname_exprs)
        } else if (length(cur_roi_cellname_indices) == 1) {
            cur_roi_cellname_means <- roi_cell_exprs[cur_roi_cellname_indices, ]
        } else {
            cur_roi_cellname_means <- rep(0, length(function_markers))
        }
        cell_function_names <- paste(cur_cn_name, function_markers, sep="-")
        names(cur_roi_cellname_means) <- cell_function_names
        roi_state_feas <- append(roi_state_feas, cur_roi_cellname_means)    
    }
    roi_fea_df[ind, ] <- roi_state_feas
}

# generate function feature names
function_fea_names <- c()
for (cur_cn_name in all_cn_lst) {
    cell_function_names <- paste(cur_cn_name, function_markers, sep="-")
    function_fea_names <- append(function_fea_names, cell_function_names)
}
colnames(roi_fea_df) <- function_fea_names
roi_fea_df <- cbind(roi_lst, roi_fea_df)
colnames(roi_fea_df)[1] <- "ROI_ID"
# save cell state features
cell_state_csv_path <- file.path(feature_dir, "CN_StateFeas.csv")
write.csv(roi_fea_df, file=cell_state_csv_path, row.names = FALSE)
