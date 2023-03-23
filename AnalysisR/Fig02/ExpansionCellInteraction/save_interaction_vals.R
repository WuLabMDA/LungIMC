library(imcRtools)
library(ggplot2)
library(tidyverse)
library(scales)
library(openxlsx)
library(stringr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
celltype_expansion_dir <- file.path(phenotype_dir, "ExpansionInteraction")
threshold_val <- 32
cell_type_interaction_name <- paste0("ExpansionInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_expansion_dir, cell_type_interaction_name)
cell_spatial_path <- file.path(cell_type_interaction_path)
load(cell_spatial_path)

# group interaction
group_interaction <- interaction_out %>% as_tibble() %>% group_by(from_label, to_label) 

## save to csv
interaction_df <- data.frame(
    roi_names = group_interaction$group_by,
    from_phenotypes = group_interaction$from_label,
    to_phenotypes = group_interaction$to_label,
    sig_vals = group_interaction$sigval
)
interaction_sigval_name <- paste0("ExpansionInteraction", threshold_val, "-Sigval.csv")
interaction_sigval_csv_path <- file.path(celltype_expansion_dir, interaction_sigval_name)
write.csv(interaction_df, file=interaction_sigval_csv_path, row.names = FALSE)