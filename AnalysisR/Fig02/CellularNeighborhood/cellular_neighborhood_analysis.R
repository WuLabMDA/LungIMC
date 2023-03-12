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

spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

set.seed(1234)
spe <- aggregateNeighbors(spe, colPairName = "neighborhood", aggregate_by = "metadata", count_by = "celltype")
cn_kmeans <- kmeans(spe$aggregatedNeighbors, centers = 10)
spe$cn_celltypes <- as.factor(cn_kmeans$cluster)

# # save the spe with cellular neighborhood information
# spe_cn_path <- file.path(data_root_dir, "SpatialAnalysis", "spe_cn10.rds")
# saveRDS(spe, spe_cn_path)

# save cn
cell_df <- data.frame(cell_id = rownames(spe@colData),
                      cell_type = spe$celltype,
                      cell_cn = spe$cn_celltypes)
cell_type_path <- file.path(phenotype_dir, "cell_cn_types.csv")
write.csv(cell_df, cell_type_path, row.names=FALSE)






