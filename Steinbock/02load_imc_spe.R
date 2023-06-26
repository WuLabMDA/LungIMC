library(imcRtools)

# set spe file path
data_root_dir <- "E:/LungIMCData/HumanSampling35/LungROIProcessing/Steinbock"
spe_file_path <- file.path(data_root_dir, "steinbock_spe.rds")

# load spe data
spe <- readRDS(spe_file_path)