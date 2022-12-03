library(imcRtools)
library(cytomapper)
library(BiocParallel)
library(tiff)


# correct all ROIs
# spillover_roi_root <- "E:/LungIMCData/HumanWholeIMC/LungROIProcessing/Denoise/SpilloverCorrection"
spillover_roi_root <- "E:/LungIMCData/TonsilIMC/TonsilROIProcessing/Denoise/SpilloverCorrection"


raw_img_dir <- file.path(spillover_roi_root, "Raw")
correct_img_dir <- file.path(spillover_roi_root, "Correct")
if (!dir.exists(correct_img_dir)) {
    dir.create(correct_img_dir)
}

# load spillover matrix
spillover_root <- "E:/LungIMCData"
spillover_matrix_path <- file.path(spillover_root, "SpilloverMatrix", "SpilloverMatrix.rds")

meta_sm <- readRDS(file = spillover_matrix_path)

# update channel names
meta_list <- c("La139", "Pr141", "Nd143", "Nd144", "Nd145", "Nd146", "Sm147",
               "Nd148", "Sm149", "Nd150", "Eu151", "Sm152", "Eu153", "Sm154",
               "Gd155", "Gd156", "Tb159", "Gd160", "Dy161", "Dy162", "Dy163",
               "Dy164", "Ho165", "Er166", "Er167", "Er168", "Tm169", "Er170",
               "Yb171", "Yb172", "Yb173", "Yb174", "Lu175", "Yb176")
channel_names <- paste0(meta_list, "Di")

file_list <- list.files(path=raw_img_dir, pattern=".tiff", all.files=TRUE, full.names=TRUE)
for (ind in 1:length(file_list)){
    start_time <- Sys.time()
    print(paste0(ind, "/", length(file_list)))
    # load
    raw_images <- loadImages(file_list[[ind]], as.is = TRUE)
    # correct
    channelNames(raw_images) <- channel_names
    correct_images <- compImage(raw_images, meta_sm, BPPARAM = SnowParam())
    # save
    roi_names <- names(correct_images)
    for (cur_name in roi_names) {
        cur_img <- as.array(correct_images[[cur_name]]) / (2^16 - 1)
        # mode(cur_img) <- "integer"
        correct_img_path <- file.path(correct_img_dir, paste0(cur_name, ".tiff"))
        writeImage(cur_img, correct_img_path, bits.per.sample=16) 
    }
    end_time <- Sys.time()
    time_dif <- difftime(end_time, start_time, units='mins')
    print(time_dif)    
}
