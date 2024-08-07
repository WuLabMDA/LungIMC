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

# detecting all cells within a given distance to the center cell (expansion)
threshold_val <- 32
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "expansion", threshold = threshold_val)
## test cell_type interaction
print(paste("Start @", format(Sys.time(), "%a %b %d %X %Y")))
interaction_out <- testInteractions(spe, group_by = "sample_id",
                                    label = "celltype", 
                                    colPairName = "expansion_interaction_graph", 
                                    method = "patch",
                                    patch_size = 3,
                                    iter = 200)
# save the results
celltype_expansion_dir <- file.path(phenotype_dir, "ExpansionInteraction")
if (!dir.exists(celltype_expansion_dir))
    dir.create(celltype_expansion_dir)
cell_type_interaction_name <- paste0("ExpansionInteractionThreshold", threshold_val, ".RData")
cell_type_interaction_path <- file.path(celltype_expansion_dir, cell_type_interaction_name)
save(spe, interaction_out, file = cell_type_interaction_path)
print(paste("Finish @", format(Sys.time(), "%a %b %d %X %Y")))