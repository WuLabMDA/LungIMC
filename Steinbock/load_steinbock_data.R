library(imcRtools)
library(cytomapper)
library(scater)
library(stringr)
library(dittoSeq)
library(RColorBrewer)
library(openxlsx)

# data_root_dir <- "E:/LungIMCData/HumanSampling41/LungROIProcessing/Steinbock"
# data_root_dir <- "E:/LungIMCData/HumanSampling35-0/LungROIProcessing/Steinbock"
# data_root_dir <- "E:/LungIMCData/HumanSampling35-1/LungROIProcessing/Steinbock"
data_root_dir <- "E:/LungIMCData/HumanSampling35-2/LungROIProcessing/Steinbock"

# Load steinbock generated data
spe <- read_steinbock(data_root_dir)

# Set the colnames of the object to generate unique identifiers per cell
colnames(spe) <- paste0(spe$sample_id, "_", spe$ObjectNumber)
# Transform signal via arcsin hyperbolic
assay(spe, "exprs") <- asinh(counts(spe)/5.0)
# Define interesting channels
rowData(spe)$use_channel <- !grepl("Ir191|NaKATPase|B2M", rownames(spe))

# add slide and roi information
spe$slide_id <- str_extract_all(spe$sample_id, ".+(?=-ROI)", simplify = TRUE)
spe$roi_id <- str_extract_all(spe$sample_id, "ROI[0-9]{3}", simplify = TRUE)

# Save generated SingleCellExperiment object
spe_file_path <- file.path(data_root_dir, "steinbock_spe.rds")
saveRDS(spe, spe_file_path)
