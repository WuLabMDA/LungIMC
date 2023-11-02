library(imcRtools)

# set spe file path
data_root_dir <- "E:/LungIMCData/HumanWholeIMC/LungROIProcessing/Steinbock"
spe_file_path <- file.path(data_root_dir, "steinbock_spe.rds")

# load spe data
spe <- readRDS(spe_file_path)