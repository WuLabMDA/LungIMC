library(imcRtools)
library(tidyverse)
library(readxl)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

## save updated data
updated_spe_name <- "lung_spe_32_cell_subtypes_final"
updated_spe_path <- file.path(phenotype_dir, paste0(updated_spe_name, ".rds"))

spe <- readRDS(updated_spe_path)
celltypes <- spe$celltype
print(paste("Total cell number: ", length(celltypes)))

img_file_path <- file.path(data_root_dir, "LungROIProcessing", "Steinbock", "images.csv")
img_data_df <- read_csv(img_file_path)
num_img <- nrow(img_data_df)
roi_cell_num <- vector("integer", num_img)

for (ind in 1:num_img) {
    img_fullname <- img_data_df$image[ind]
    img_name <- substr(img_fullname, 1, nchar(img_fullname) - 5)
    roi_object <- spe[, spe$sample_id == img_name]
    roi_cell_num[ind] <- length(roi_object$celltype)
}