library(SpatialExperiment)
library(imcRtools)
library(BiocParallel)
library(cytomapper)
library(tidyverse)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

## save to csv
cell_id_lst <- rownames(colData(spe))
celltype_lst <- spe$celltype
cell_df <- data.frame(
    cell_id = cell_id_lst,
    celltype = celltype_lst
)
celltype_csv_path <- file.path(phenotype_dir, "cell_phenotypes.csv")
write.csv(cell_df, file=celltype_csv_path, row.names = FALSE)