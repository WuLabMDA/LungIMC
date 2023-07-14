library(imcRtools)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

## load raw data
raw_spe_name <- "steinbock_spe-all_celltypes_Renamed"
raw_spe_path <- file.path(phenotype_dir, paste0(raw_spe_name, ".rds"))
spe <- readRDS(raw_spe_path)

### Update cell type information
renamed_celltypes <- spe@metadata$celltypesRenamed
renamed_celltypes[renamed_celltypes == "T-reg-Cell"] <- "T-Reg-Cell"
spe$celltype <- renamed_celltypes

## save updated data
updated_spe_name <- "lung_spe_15_celltypes_final"
updated_spe_path <- file.path(phenotype_dir, paste0(updated_spe_name, ".rds"))
saveRDS(spe, file = updated_spe_path)