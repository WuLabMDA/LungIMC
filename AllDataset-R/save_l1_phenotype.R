library(SpatialExperiment)

set.seed(1234)
data_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockAll"
phenotype_path <- file.path(data_root_dir, "som_lineage_spe.rds")
spe <- readRDS(phenotype_path)

# collect cell id and phenotype information
cell_ids <- colnames(spe)
cell_phenotypes <- spe$som_clusters

# save the cell ids with matching prototypes
cell_id_phenotype_df <- data.frame(CellID = cell_ids, CellPhenotype = cell_phenotypes)
cell_id_phenotype_csv_path <- file.path(data_root_dir, "cell_id_phenotype.csv")
write.csv(cell_id_phenotype_df, cell_id_phenotype_csv_path, row.names = FALSE)