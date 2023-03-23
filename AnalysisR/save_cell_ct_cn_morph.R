library(SpatialExperiment)
library(imcRtools)
library(cytomapper)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(viridis)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(readxl)
library(knitr)

## load all interactions
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

# load the spe
threshold_val <- 50
num_cn <- 10
spe_cn_name <- paste0("Delaunay", threshold_val, "-CN", num_cn, ".rds")
spe_cn_path <- file.path(phenotype_dir, spe_cn_name)
spe <- readRDS(spe_cn_path)

# save cell information
cell_df <- data.frame(cell_id = rownames(spe@colData),
                      cell_type = spe$celltype,
                      cell_cn = spe$cn_celltypes,
                      cell_area = spe$area,
                      cell_maj_ax_len = spe$major_axis_length,
                      cell_min_ax_len = spe$minor_axis_length,
                      cell_eccentricity = spe$eccentricity)
cell_type_path <- file.path(phenotype_dir, "cell_ct_cn_morphs.csv")
write.csv(cell_df, cell_type_path, row.names=FALSE)






