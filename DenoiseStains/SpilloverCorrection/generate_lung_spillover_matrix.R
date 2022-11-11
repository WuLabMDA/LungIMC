library(imcRtools)
library(CATALYST)
library(pheatmap)


spillover_root <- "E:/LungIMCData/LungROIProcessing/Denoise/SpilloverCorrection"

# load spots data
dset <- "LungIMC"
spots_path <- file.path(spillover_root, "SpilloverMatrix", dset)
sce <- readSCEfromTXT(spots_path) 
# plotSpotHeatmap(sce)

# pixel binning - optional
bin_size = 10
sce <- binAcrossPixels(sce, bin_size = bin_size)

# asinh-transform the data using a cofactor of 5
assay(sce, "exprs") <- asinh(counts(sce)/5)

# Filtering incorrectly assigned pixels
bc_key <- as.numeric(unique(sce$sample_mass))
bc_key <- bc_key[order(bc_key)]
sce <- assignPrelim(sce, bc_key = bc_key)
sce <- estCutoffs(sce)
sce <- applyCutoffs(sce)


# Exclude incorrect assigned pixels
sce <- filterPixels(sce, minevents = 40, correct_pixels = TRUE)

# Compute spillover matrix
sce <- computeSpillmat(sce)
meta_list <- c("La139", "Pr141", "Nd143", "Nd144", "Nd145", "Nd146", "Sm147",
               "Nd148", "Sm149", "Nd150", "Eu151", "Sm152", "Eu153", "Sm154",
               "Gd155", "Gd156", "Tb159", "Gd160", "Dy161", "Dy162", "Dy163",
               "Dy164", "Ho165", "Er166", "Er167", "Er168", "Tm169", "Er170",
               "Yb171", "Yb172", "Yb173", "Yb174", "Lu175", "Yb176")
channel_names <- paste0(meta_list, "Di")
sm <- metadata(sce)$spillover_matrix
meta_sm <- adaptSpillmat(sm, out_chs = channel_names)

spillover_matrix_path <- file.path(spillover_root, "SpilloverMatrix", "SpilloverMatrix.rds")
saveRDS(meta_sm, file = spillover_matrix_path)

# rowData(sce)$is_bc <- TRUE # to make sure all antibodies show up
# plotSpillmat(sce, sm=meta_sm)
