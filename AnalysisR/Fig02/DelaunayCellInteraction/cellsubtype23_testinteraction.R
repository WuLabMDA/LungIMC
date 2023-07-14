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
spe_cellsubtype_name <-"lung_spe_33_cell_subtypes_final"
spe_cellsubtype_path <- file.path(phenotype_dir, paste0(spe_cellsubtype_name, ".rds"))
spe <- readRDS(spe_cellsubtype_path)

# merge cell subtypes
spe$cellsubtype[spe$cellsubtype=="Ki67+ Epithelial"] <- "Epithelial"
spe$cellsubtype[spe$cellsubtype=="Ki67+ PDL1+ Epithelial"] <- "Epithelial"
spe$cellsubtype[spe$cellsubtype=="CD73+ Epithelial"] <- "Epithelial"
spe$cellsubtype[spe$cellsubtype=="Other Epithelial"] <- "Epithelial"
spe$cellsubtype[spe$cellsubtype=="Ki67+ B-Cells"] <- "B-Cells"
spe$cellsubtype[spe$cellsubtype=="Ki67- B-Cells"] <- "B-Cells"
spe$cellsubtype[spe$cellsubtype=="Ki67+ NK Cells"] <- "NK-Cells"
spe$cellsubtype[spe$cellsubtype=="Ki67- NK Cells"] <- "NK-Cells"
spe$cellsubtype[spe$cellsubtype=="Ki67+ Dendritic Cells"] <- "Dendritic-Cell"
spe$cellsubtype[spe$cellsubtype=="HLADR+ Dendritic Cells"] <- "Dendritic-Cell"
spe$cellsubtype[spe$cellsubtype=="Other Dendritic Cells"] <- "Dendritic-Cell"
spe$cellsubtype[spe$cellsubtype=="Ki67+ Treg-Cells"] <- "Treg-Cells"
spe$cellsubtype[spe$cellsubtype=="Ki67- Treg-Cells"] <- "Treg-Cells"
spe$cellsubtype[spe$cellsubtype=="CD163+ Ki67+ Macrophage"] <- "CD163+ Macrophage"
spe$cellsubtype[spe$cellsubtype=="CD163- PDL1+ Macrophage"] <- "CD163- Macrophage"

# detecting all cells within a given distance to the center cell (expansion)
threshold_val <- 50
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "delaunay", max_dist = threshold_val)

## test cell_type interaction
print(paste("Start @", format(Sys.time(), "%a %b %d %X %Y")))
interaction_out <- testInteractions(spe, group_by = "sample_id",
                                    label = "cellsubtype", 
                                    colPairName = "delaunay_interaction_graph", 
                                    method = "patch",
                                    patch_size = 3,
                                    iter = 200)
# save the results
celltype_delaunay_dir <- file.path(phenotype_dir, "DelaunayInteraction")
if (!dir.exists(celltype_delaunay_dir))
    dir.create(celltype_delaunay_dir)

cell_subtype_interaction_name <- paste0("Subtype23DelaunayInteractionThreshold", threshold_val, ".RData")
cell_subtype_interaction_path <- file.path(celltype_delaunay_dir, cell_subtype_interaction_name)
save(spe, interaction_out, file = cell_subtype_interaction_path)
print(paste("Finish @", format(Sys.time(), "%a %b %d %X %Y")))