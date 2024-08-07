library(imcRtools)
library(readxl)
library(stringr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(gridExtra)
library(corrplot)
library(RColorBrewer)
library(NCmisc)
library(dittoSeq)
library(viridis)
library(hrbrthemes)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

spe_celltype_name <-"lung_spe_32_cell_subtypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

metadata_dir <- file.path(data_root_dir, "Metadata")
# load ROI information
roi_info_name <- "ROI_Info_Aggregation.csv"
roi_info_path <- file.path(metadata_dir, roi_info_name)
roi_df <- read_csv(roi_info_path)
# load lesion information
slide_info_name <- "Lesion_Info_Aggregation.csv"
slide_info_path <- file.path(metadata_dir, slide_info_name)
slide_info_df <- read_csv(slide_info_path)
# load patient information
patient_info_name <- "Patient_Info.xlsx"
patient_info_path <- file.path(metadata_dir, patient_info_name)
patient_df <- read_excel(patient_info_path)

## filtering ROIs
roi_lst <- roi_df$ROI_ID
interested_roi_vec <- c()
for (ind in 1:nrow(roi_df)) {
    roi_location <- roi_df$ROI_Location[ind]
    if (roi_location %in% c("Normal", "Tumor")) 
        interested_roi_vec <- append(interested_roi_vec, roi_df$ROI_ID[ind])
}

## obtain cells inside tumor lesions
lesion_indices = c()
cell_id_lst <- rownames(colData(spe))
for (cur_roi in interested_roi_vec) 
    lesion_indices <- append(lesion_indices, which(startsWith(cell_id_lst, cur_roi)))
cell_id_lst <- cell_id_lst[lesion_indices]
celltype_lst <- spe$cellsubtype
celltype_lst <- celltype_lst[lesion_indices]


subtype_proportion_dir <- file.path(data_root_dir, "NatureFigures", "Fig02", "ProportionDensity", "SubtypeProportion")
if (!file.exists(subtype_proportion_dir))
    dir.create(subtype_proportion_dir, recursive = TRUE)

m2_cell_types <- c("CD163+ Macrophages", "Ki67+ Macrophages")
m1_cell_types <- c("CD163- Macrophages")
# collect information
ratio_lst <- c()
stage_lst <- c()

# iterate by slide
stage_slide_lst <- slide_info_df$Slide_ID
for (ir in 1:length(stage_slide_lst)) {
    # iterate by slide
    cur_slide <- stage_slide_lst[ir]
    if (cur_slide == "2571-1D")
        next
    slide_cell_indices <- which(startsWith(cell_id_lst, cur_slide))
    slide_celltypes <- celltype_lst[slide_cell_indices]
    m2_cell_num <- sum(slide_celltypes %in% m2_cell_types)
    m1_cell_num <- sum(slide_celltypes %in% m1_cell_types)
    slide_cell_ratio <- (m2_cell_num + 0.01) / (m1_cell_num + 0.01) 
    # gather information
    ratio_lst <- append(ratio_lst, slide_cell_ratio)
    slide_index <- which(slide_info_df$Slide_ID == cur_slide)
    stage_lst <- append(stage_lst, slide_info_df$Slide_Diag[slide_index])
}

cell_ratio_df <- data.frame(Ratio=ratio_lst, Stage=stage_lst)
long_ratio_df <- gather(cell_ratio_df, key = "variable", value = "value", -Ratio)

cell_ratio_df$Stage <- factor(cell_ratio_df$Stage, levels = c("Normal", "AAH", "AIS", "MIA", "ADC"))
ggplot(cell_ratio_df, aes(x = Stage, y = Ratio, colour = Stage)) + 
    geom_boxplot() +
    geom_jitter(width = 0.25) + 
    theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle(paste0("M21-Ratio", " - Proportion"))

subtype_proportion_plot_path <- file.path(subtype_proportion_dir, paste0("M21-Ratio", ".pdf"))
ggsave(filename = subtype_proportion_plot_path, device='pdf', width=5, height=8, dpi=300)
while (!is.null(dev.list()))
    dev.off()    

subtype_proportion_pfile <- file(file.path(subtype_proportion_dir, paste0("M21-Ratio", "-pval.txt")), "w")
# Welch Two Sample t-test p-value
normal_aah_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Normal"], long_ratio_df$Ratio[long_ratio_df$value=="AAH"])
writeLines(paste0("Normal-AAH:", normal_aah_ttest$p.value), subtype_proportion_pfile)
normal_ais_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Normal"], long_ratio_df$Ratio[long_ratio_df$value=="AIS"])
writeLines(paste0("Normal-AIS:", normal_ais_ttest$p.value), subtype_proportion_pfile)
normal_mia_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Normal"], long_ratio_df$Ratio[long_ratio_df$value=="MIA"])
writeLines(paste0("Normal-MIA:", normal_mia_ttest$p.value), subtype_proportion_pfile)
normal_adc_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="Normal"], long_ratio_df$Ratio[long_ratio_df$value=="ADC"])
writeLines(paste0("Normal-ADC:", normal_adc_ttest$p.value), subtype_proportion_pfile)
aah_ais_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AAH"], long_ratio_df$Ratio[long_ratio_df$value=="AIS"])
writeLines(paste0("AAH-AIS:", aah_ais_ttest$p.value), subtype_proportion_pfile)
aah_mia_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AAH"], long_ratio_df$Ratio[long_ratio_df$value=="MIA"])
writeLines(paste0("AAH-MIA:", aah_mia_ttest$p.value), subtype_proportion_pfile)
aah_adc_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AAH"], long_ratio_df$Ratio[long_ratio_df$value=="ADC"])
writeLines(paste0("AAH-ADC:", aah_adc_ttest$p.value), subtype_proportion_pfile)
ais_mia_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AIS"], long_ratio_df$Ratio[long_ratio_df$value=="MIA"])
writeLines(paste0("AIS-MIA:", ais_mia_ttest$p.value), subtype_proportion_pfile)
ais_adc_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="AIS"], long_ratio_df$Ratio[long_ratio_df$value=="ADC"])
writeLines(paste0("AIS-ADC:", ais_adc_ttest$p.value), subtype_proportion_pfile)
mia_adc_ttest <- t.test(long_ratio_df$Ratio[long_ratio_df$value=="MIA"], long_ratio_df$Ratio[long_ratio_df$value=="ADC"])
writeLines(paste0("MIA-ADC:", mia_adc_ttest$p.value), subtype_proportion_pfile)
close(subtype_proportion_pfile)
