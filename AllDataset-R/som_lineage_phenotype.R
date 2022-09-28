library(SpatialExperiment)
library(CATALYST)
library(scran)
library(bluster)
library(BiocParallel)
library(kohonen)
library(ConsensusClusterPlus)
library(dittoSeq)
library(viridis)
set.seed(1234)


data_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockAll"

raw_spe_path = file.path(data_root_dir, "all_steinbock_spe.rds")
spe <- readRDS(raw_spe_path)

# Define lineage marker list
lineage_strings <- "CD45|CD3e|CD4|FoxP3|CD8a|CD19|CD94|CD11b|CD11c|CD14|MPO|CD68|CD33|HLADR|CD45RO|CK|CD31|aSMA"
rowData(spe)$lineage_markers <- grepl(lineage_strings, rownames(spe))

# Run FlowSOM and ConsensusClusterPlus clustering
spe <- cluster(spe, features = rownames(spe)[rowData(spe)$lineage_markers], maxK = 30, seed = 1234)

# # Assess cluster stability
# delta_area(spe)

# save som clustered results
spe$som_clusters <- cluster_ids(spe, "meta10")
phenotype_path <- file.path(data_root_dir, "som_raw_spe.rds")
saveRDS(spe, phenotype_path)