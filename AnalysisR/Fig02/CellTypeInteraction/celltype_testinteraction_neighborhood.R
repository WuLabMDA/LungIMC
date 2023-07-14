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

## load spe object
spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

## test cell_type interaction
print(paste("Start @", format(Sys.time(), "%a %b %d %X %Y")))
interaction_out <- testInteractions(spe, group_by = "sample_id",
                                    label = "celltype", 
                                    colPairName = "neighborhood", 
                                    method = "patch",
                                    patch_size = 3,
                                    iter = 200)
cell_type_spatial_path <- file.path(phenotype_dir, paste0("InteractionsTestIter200Neighborhood", ".RData"))
save(spe, interaction_out, file = cell_type_spatial_path)
print(paste("Finish @", format(Sys.time(), "%a %b %d %X %Y")))