library(SpatialExperiment)


data_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockAll"
all_phenotype_path <- file.path(data_root_dir, "som_lineage_spe.rds")
all_spe <- readRDS(all_phenotype_path)