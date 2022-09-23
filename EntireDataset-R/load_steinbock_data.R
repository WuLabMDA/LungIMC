library(imcRtools)
library(cytomapper)
library(scater)
library(stringr)
library(dittoSeq)
library(RColorBrewer)
library(openxlsx)

data_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockAll"
# Load steinbock generated data
spe <- read_steinbock(data_root_dir)

# Set the colnames of the object to generate unique identifiers per cell
colnames(spe) <- paste0(spe$sample_id, "_", spe$ObjectNumber)
# Transform signal via arcsin hyperbolic
assay(spe, "exprs") <- asinh(counts(spe)/5.0)
# Define interesting channels
rowData(spe)$use_channel <- !grepl("Ir191|NaKATPase|B2M", rownames(spe))
# Embed expression via UMAP
spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs")

# add slide and roi information
spe$slide_id <- str_extract_all(spe$sample_id, ".+(?=-ROI)", simplify = TRUE)
spe$roi_id <- str_extract_all(spe$sample_id, "ROI[0-9]{3}", simplify = TRUE)

# add location & diagnosis information
roi_info <- read.xlsx(file.path(data_root_dir, "StudyROI_Info.xlsx"))
spe$roi_location <- roi_info$ROI_Location[match(spe$sample_id, roi_info$ROI_ID)]
spe$roi_diag <- roi_info$ROI_Diag[match(spe$sample_id, roi_info$ROI_ID)]

# add slide-level patient and batch information
slide_info <- read.xlsx(file.path(data_root_dir, "StudySlide_Info.xlsx"))
spe$slide_diag <- slide_info$Slide_Diag[match(spe$slide_id, slide_info$Slide_ID)]
spe$patient_id <- slide_info$Patient_ID[match(spe$slide_id, slide_info$Slide_ID)]
spe$stain_id <- slide_info$Staining_ID[match(spe$slide_id, slide_info$Slide_ID)]

# Define color schemes
color_vectors <- list()
color_vectors$ROILoc <- setNames(brewer.pal(length(unique(spe$roi_location)), name = "Set1"), unique(spe$roi_location))
color_vectors$ROIDiag <- setNames(brewer.pal(length(unique(spe$roi_diag)), name = "Set1"), unique(spe$roi_diag))
color_vectors$SlideDiag <- setNames(brewer.pal(length(unique(spe$slide_diag)), name = "Set1"), unique(spe$slide_diag))
color_vectors$PatientID <- setNames(colorRampPalette(brewer.pal(4, "PuOr"))(length(unique(spe$patient_id))), unique(spe$patient_id))
color_vectors$Stain <- setNames(colorRampPalette(brewer.pal(4, "PuOr"))(length(unique(spe$stain_id))), unique(spe$stain_id))
metadata(spe)$color_vectors <- color_vectors

# Save generated SingleCellExperiment object
raw_spe_path = file.path(data_root_dir, "all_steinbock_spe.rds")
saveRDS(spe, raw_spe_path)