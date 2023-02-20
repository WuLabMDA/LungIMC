library(imcRtools)
library(cytomapper)
library(stringr)
library(readxl)
library(RColorBrewer)

# set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
## create spatial analysis directory
spatial_dir <- file.path(data_root_dir, "SpatialAnalysis")
if (!dir.exists(spatial_dir))
    dir.create(spatial_dir)
## load phenotype information
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

# add meta information
metadata_dir <- file.path(data_root_dir, "Metadata")
## Add ROI information
roi_info_name <- "ROI_Info"
roi_info_path <- file.path(metadata_dir, paste0(roi_info_name, ".xlsx"))
roi_df <- read_excel(roi_info_path)
spe$roi_loc <- roi_df$ROI_Location[match(spe$roi_id, roi_df$ROI_ID)]
spe$roi_diag <- roi_df$ROI_Diag[match(spe$roi_id, roi_df$ROI_ID)]
## Add lesion information
slide_info_name <- "Lesion_Info"
slide_info_path <- file.path(metadata_dir, paste0(slide_info_name, ".xlsx"))
slide_df <- read_excel(slide_info_path)
spe$slide_diag <- slide_df$Slide_Diag[match(spe$slide_id, slide_df$Slide_ID)]
spe$slide_pat <- slide_df$Patient_ID[match(spe$slide_id, slide_df$Slide_ID)]
## Add patient information
patient_match_regex <- ".+(?=-[[0-9]{1,}|[a-zA-Z]{1,}][[0-9]|[a-zA-Z]]*$)"
spe$patient_id <- str_extract_all(spe$slide_id, patient_match_regex, simplify = TRUE)
patient_info_name <- "Patient_Info"
patient_info_path <- file.path(metadata_dir, paste0(patient_info_name, ".xlsx"))
patient_df <- read_excel(patient_info_path)
spe$pat_gender <- patient_df$Gender[match(spe$patient_id, patient_df$PatientID)]
spe$pat_race <- patient_df$Race[match(spe$patient_id, patient_df$PatientID)]
patient_df$RecurrentStatus[patient_df$RecurrentStatus == "RECURRENT"] <- "Recur"
patient_df$RecurrentStatus[patient_df$RecurrentStatus == "NO RECURRENT"] <- "Non-Recur"
spe$pat_recurrence <- patient_df$RecurrentStatus[match(spe$patient_id, patient_df$PatientID)]
## median_age <- median(patient_df$Age)
median_age <- 71.3
patient_df$Age[patient_df$Age <= median_age] <- "<=71"
patient_df$Age[patient_df$Age > median_age] <- ">71"
spe$pat_age <- patient_df$Age[match(spe$patient_id, patient_df$PatientID)]
spe$pat_smoke <- patient_df$SmokeStatus[match(spe$patient_id, patient_df$PatientID)]

## save updated spe
updated_spe_name <- "lung_imc_spe"
updated_spe_path <- file.path(spatial_dir, paste0(updated_spe_name, ".rds"))
saveRDS(spe, file = updated_spe_path)

steinbock_dir <- file.path(data_root_dir, "LungROIProcessing", "Steinbock")
# load images
imc_img_dir <- file.path(steinbock_dir, "img")
images <- loadImages(imc_img_dir)
channelNames(images) <- rownames(spe)
# load masks
imc_mask_dir <- file.path(steinbock_dir, "masks_deepcell")
masks <- loadImages(imc_mask_dir, as.is = TRUE)
# all.equal(names(images), names(masks))

img_slide_ids <- str_extract_all(names(masks), ".+(?=-ROI)", simplify = TRUE)
img_pat_id <- str_extract_all(img_slide_ids, patient_match_regex, simplify = TRUE)
img_stages <- roi_df$ROI_Diag[match(names(masks), roi_df$ROI_ID)]
img_info_df <- DataFrame(sample_id = names(masks), patient_id = img_pat_id, stage = img_stages)
mcols(images) <- img_info_df
mcols(masks) <- img_info_df
# save images
lung_imc_img_path <- file.path(spatial_dir, "lung_imc_images.rds")
saveRDS(images, lung_imc_img_path)
# save masks
lung_imc_mask_path <- file.path(spatial_dir, "lung_imc_masks.rds")
saveRDS(masks, lung_imc_mask_path)


