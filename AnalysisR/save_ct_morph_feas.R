library(SpatialExperiment)
library(imcRtools)
library(BiocParallel)
library(cytomapper)
library(tidyverse)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
feature_dir <- file.path(data_root_dir, "FeatureAnalysis")

spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

# extract cell information
cell_id_lst <- rownames(colData(spe))
cell_phenotype_lst <- spe$celltype

all_cell_lst <- c("Epithelial-Cell", "Endothelial-Cell", "Fibroblast", "CD4-T-Cell", "CD8-T-Cell", 
                  "T-Reg-Cell", "B-Cell", "Macrophage", "Monocyte", "Dendritic-Cell", 
                  "Neutrophil", "MDSC", "NK-Cell", "Proliferating-Cell", "Undefined")

# obtain morphology features
morph_fea_names <- c("Area", "MajorAxisLength", "MinorAxisLength", "Eccentricity")
morph_df <- data.frame(cell_area = spe$area,
                       cell_maj_ax_len = spe$major_axis_length,
                       cell_min_ax_len = spe$minor_axis_length,
                       cell_eccentricity = spe$eccentricity)
morph_mat <- data.matrix(morph_df)
rownames(morph_mat) <- rownames(spe@colData)
colnames(morph_mat) <- morph_fea_names


roi_lst <- unique(spe$sample_id)
morph_fea_num <- length(all_cell_lst) * length(morph_fea_names)
roi_fea_df <- data.frame(matrix(ncol = morph_fea_num, nrow = length(roi_lst)))
# traverse ROI one-by-one
for (ind in 1:length(roi_lst)) {
    # extract cell index in each ROI
    cur_roi <- roi_lst[ind]
    roi_cell_indices <- startsWith(cell_id_lst, cur_roi)
    roi_cell_phenotypes <- cell_phenotype_lst[roi_cell_indices]
    roi_cell_morphs <- morph_mat[roi_cell_indices, ]
    # aggreate functional marker states
    roi_morph_feas <- c()
    for (cur_cell_name in all_cell_lst) {
        cur_roi_cellname_indices <- which(roi_cell_phenotypes == cur_cell_name)
        if (length(cur_roi_cellname_indices) > 1) {
            cur_roi_cellname_morphs <- roi_cell_morphs[cur_roi_cellname_indices, ]
            cur_roi_cellname_means <- colMeans(cur_roi_cellname_morphs)
        } else if (length(cur_roi_cellname_indices) == 1) {
            cur_roi_cellname_means <- roi_cell_morphs[cur_roi_cellname_indices, ]
        } else {
            cur_roi_cellname_means <- rep(0, length(morph_fea_names))
        }
        cell_morph_names <- paste(cur_cell_name, morph_fea_names, sep="-")
        names(cur_roi_cellname_means) <- cell_morph_names
        roi_morph_feas <- append(roi_morph_feas, cur_roi_cellname_means)    
    }
    roi_fea_df[ind, ] <- roi_morph_feas
}


# generate morph feature names
all_morph_fea_names <- c()
for (cur_cell_name in all_cell_lst) {
    cell_morph_names <- paste(cur_cell_name, morph_fea_names, sep="-")
    all_morph_fea_names <- append(all_morph_fea_names, cell_morph_names)
}
colnames(roi_fea_df) <- all_morph_fea_names
roi_fea_df <- cbind(roi_lst, roi_fea_df)
colnames(roi_fea_df)[1] <- "ROI_ID"
# save cell state features
cell_morph_csv_path <- file.path(feature_dir, "CT_MorphFeas.csv")
write.csv(roi_fea_df, file=cell_morph_csv_path, row.names = FALSE)
