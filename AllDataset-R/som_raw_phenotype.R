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


normal_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockAll"

raw_spe_path = file.path(normal_root_dir, "all_steinbock_spe.rds")
spe <- readRDS(raw_spe_path)

# Run FlowSOM and ConsensusClusterPlus clustering
spe <- cluster(spe, features = rownames(spe)[rowData(spe)$use_channel], maxK = 20, seed = 1234)

# # Assess cluster stability
# delta_area(spe)

# save som clustered results
spe$som_clusters <- cluster_ids(spe, "meta12")
phenotype_path <- file.path(normal_root_dir, "som_raw_spe.rds")
saveRDS(spe, phenotype_path)