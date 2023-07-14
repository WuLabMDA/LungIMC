library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
celltype_delaunay_dir <- file.path(phenotype_dir, "DelaunayInteraction")
# threshold_val <- 30
# threshold_val <- 40
threshold_val <- 50
cell_type_interaction_name <- paste0("DelaunayInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_delaunay_dir, cell_type_interaction_name)
load(cell_type_interaction_path)

# replace NA to -1, both unavailable 0
cell_type_num <- 15
interaction_num <- cell_type_num * cell_type_num
num_roi <- interaction_out@nrows / interaction_num
for (nn in 1:num_roi) {
    start_ind <- (nn - 1) * interaction_num + 1
    end_ind <- nn * interaction_num
    cur_sigvals <- interaction_out$sigval[start_ind:end_ind]
    cell_exists <- rep(0, cell_type_num)
    for (mm in 1:cell_type_num) {
        cur_cell_vals <- cur_sigvals[(mm-1)*cell_type_num+1:mm*cell_type_num]
        if (all(is.na(cur_cell_vals)))
            cell_exists[mm] <- 1
    }
    cur_sigvals <- replace(cur_sigvals, is.na(cur_sigvals), -1)
    missing_inds <- which(cell_exists == 1)
    if (length(missing_inds) > 0) {
        ind_combs <- crossing(missing_inds, missing_inds)
        for (jj in 1:dim(ind_combs)[1]) {
            ele_ind <- (ind_combs[[jj,1]]-1)* cell_type_num + ind_combs[[jj,2]]
            cur_sigvals[ele_ind] <- 0
        }
    }
    interaction_out$sigval[start_ind:end_ind] <- cur_sigvals
}

# group interaction
group_interaction <- interaction_out %>% as_tibble() %>% group_by(from_label, to_label) 

## save to csv
interaction_df <- data.frame(
    roi_names = group_interaction$group_by,
    from_phenotypes = group_interaction$from_label,
    to_phenotypes = group_interaction$to_label,
    sig_vals = group_interaction$sigval
)

interaction_sigval_name <- paste0("DelaunayInteraction", threshold_val, "-Sigval.csv")
interaction_sigval_csv_path <- file.path(celltype_delaunay_dir, interaction_sigval_name)
write.csv(interaction_df, file=interaction_sigval_csv_path, row.names = FALSE)