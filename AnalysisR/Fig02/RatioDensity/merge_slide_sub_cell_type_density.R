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

spe_celltype_name <-"lung_spe_33_cell_subtypes_final"
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


subtype_density_dir <- file.path(data_root_dir, "NatureFigures", "Fig02", "ProportionDensity", "SubtypeDensity")
if (!file.exists(subtype_density_dir))
    dir.create(subtype_density_dir, recursive = TRUE)

# list all cell subtypes
all_cell_lst <- c("Ki67+ Epithelial", "Ki67+ PDL1+ Epithelial", "CD73+ Epithelial", "Other Epithelial",
                  "Ki67+ B-Cells", "Ki67- B-Cells", "Neutrophil", "Ki67+ NK Cells", "Ki67- NK Cells", 
                  "Ki67+ Dendritic Cells", "HLADR+ Dendritic Cells", "Other Dendritic Cells", "Endothelial-Cell",
                  "Cytotoxic CD8 T-Cells", "Memory CD8 T-Cells", "Exhausted CD8 T-Cells", "Ki67+ CD8 T-Cells", "Naive CD8 T-Cells",
                  "Memory CD4 T-Cells", "Exhausted CD4 T-Cells", "Ki67+ CD4 T-Cells", "Naive CD4 T-Cells", 
                  "Ki67+ Treg-Cells", "Ki67- Treg-Cells", "Proliferating-Cell", "CD163+ Macrophage",
                  "CD163+ Ki67+ Macrophage", "CD163- Macrophage", "CD163- PDL1+ Macrophage", "Monocyte",
                  "MDSC", "Fibroblast", "Undefined")

for (cell_type in all_cell_lst) {
    # collect information
    density_lst <- c()
    stage_lst <- c()
    
    # iterate by slide
    stage_slide_lst <- slide_info_df$Slide_ID
    for (ir in 1:length(stage_slide_lst)) {
        # iterate by slide
        cur_slide <- stage_slide_lst[ir]
        if (cur_slide == "2571-1D")
            next
        slide_cell_indices <- which(startsWith(cell_id_lst, cur_slide))
        slide_cell_rois <- unlist(unique(lapply(strsplit(cell_id_lst[slide_cell_indices], "_"), `[`, 1)))
        slide_roi_df <- roi_df[roi_df$ROI_ID %in% slide_cell_rois,]
        slide_roi_area <- sum(slide_roi_df$Area)
        slide_celltypes <- celltype_lst[slide_cell_indices]
        slide_cell_density <- sum(slide_celltypes == cell_type) / slide_roi_area
        # gather information
        density_lst <- append(density_lst, slide_cell_density)
        slide_index <- which(slide_info_df$Slide_ID == cur_slide)
        stage_lst <- append(stage_lst, slide_info_df$Slide_Diag[slide_index])    
    }   
    
    cell_density_df <- data.frame(Density=density_lst, Stage=stage_lst)
    long_density_df <- gather(cell_density_df, key = "variable", value = "value", -Density)    
    cell_density_df$Stage <- factor(cell_density_df$Stage, levels = c("Normal", "AAH", "AIS", "MIA", "ADC"))
    ggplot(cell_density_df, aes(x = Stage, y = Density, colour = Stage)) + 
        geom_boxplot() +
        geom_jitter(width = 0.25) + 
        theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        ggtitle(paste0(cell_type, " - Density"))
    
    subtype_density_plot_path <- file.path(subtype_density_dir, paste0(cell_type, ".pdf"))
    ggsave(filename = subtype_density_plot_path, device='pdf', width=5, height=8, dpi=300)
    while (!is.null(dev.list()))
        dev.off()        
    
    subtype_density_pfile <- file(file.path(subtype_density_dir, paste0(cell_type, "-pval.txt")), "w")
    normal_aah_ttest <- t.test(long_density_df$Density[long_density_df$value=="Normal"], long_density_df$Density[long_density_df$value=="AAH"])
    writeLines(paste0("Normal-AAH:", normal_aah_ttest$p.value), subtype_density_pfile)
    normal_ais_ttest <- t.test(long_density_df$Density[long_density_df$value=="Normal"], long_density_df$Density[long_density_df$value=="AIS"])
    writeLines(paste0("Normal-AIS:", normal_ais_ttest$p.value), subtype_density_pfile)
    normal_mia_ttest <- t.test(long_density_df$Density[long_density_df$value=="Normal"], long_density_df$Density[long_density_df$value=="MIA"])
    writeLines(paste0("Normal-MIA:", normal_mia_ttest$p.value), subtype_density_pfile)
    normal_adc_ttest <- t.test(long_density_df$Density[long_density_df$value=="Normal"], long_density_df$Density[long_density_df$value=="ADC"])
    writeLines(paste0("Normal-ADC:", normal_adc_ttest$p.value), subtype_density_pfile) 
    aah_ais_ttest <- t.test(long_density_df$Density[long_density_df$value=="AIS"], long_density_df$Density[long_density_df$value=="AAH"])
    writeLines(paste0("AAH-AIS:", aah_ais_ttest$p.value), subtype_density_pfile)
    aah_mia_ttest <- t.test(long_density_df$Density[long_density_df$value=="AAH"], long_density_df$Density[long_density_df$value=="MIA"])
    writeLines(paste0("AAH-MIA:", aah_mia_ttest$p.value), subtype_density_pfile)  
    aah_adc_ttest <- t.test(long_density_df$Density[long_density_df$value=="AAH"], long_density_df$Density[long_density_df$value=="ADC"])
    writeLines(paste0("AAH-ADC:", aah_adc_ttest$p.value), subtype_density_pfile)    
    ais_mia_ttest <- t.test(long_density_df$Density[long_density_df$value=="AIS"], long_density_df$Density[long_density_df$value=="MIA"])
    writeLines(paste0("AIS-MIA:", ais_mia_ttest$p.value), subtype_density_pfile)     
    ais_adc_ttest <- t.test(long_density_df$Density[long_density_df$value=="AIS"], long_density_df$Density[long_density_df$value=="ADC"])
    writeLines(paste0("AIS-ADC:", ais_adc_ttest$p.value), subtype_density_pfile)       
    mia_adc_ttest <- t.test(long_density_df$Density[long_density_df$value=="MIA"], long_density_df$Density[long_density_df$value=="ADC"])
    writeLines(paste0("MIA-ADC:", mia_adc_ttest$p.value), subtype_density_pfile)
    close(subtype_density_pfile)
}