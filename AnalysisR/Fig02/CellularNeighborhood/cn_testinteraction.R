library(imcRtools)
library(SpatialExperiment)
library(stringr)
library(RColorBrewer)
library(dittoSeq)
library(ggplot2)
library(viridis)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
celltype_delaunay_dir <- file.path(phenotype_dir, "DelaunayInteraction")

# load the spe
spe_cn_path <- file.path(celltype_delaunay_dir, "Delaunay50-CN8.rds")
spe <- readRDS(spe_cn_path)

## test cell_type interaction
print(paste("Start @", format(Sys.time(), "%a %b %d %X %Y")))
interaction_out <- testInteractions(spe, group_by = "sample_id",
                                    label = "cn_celltypes", 
                                    colPairName = "delaunay_interaction_graph", 
                                    method = "patch",
                                    patch_size = 3,
                                    iter = 200)

cn_interaction_path <- file.path(celltype_delaunay_dir, "Delaunay50-CN8-Interaction.RData")
save(spe, interaction_out, file = cn_interaction_path)
print(paste("Finish @", format(Sys.time(), "%a %b %d %X %Y")))