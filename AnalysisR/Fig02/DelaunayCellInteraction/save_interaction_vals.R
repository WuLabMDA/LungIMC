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
threshold_val <- 50
cell_type_interaction_name <- paste0("DelaunayInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_delaunay_dir, cell_type_interaction_name)
load(cell_type_interaction_path)

# replace NA to -1
interaction_out$sigval <- replace(interaction_out$sigval, is.na(interaction_out$sigval), -1)

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