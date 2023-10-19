library(imcRtools)
library(readxl)
library(ggplot2)
library(ggbeeswarm)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")

spe_celltype_name <-"lung_spe_15_celltypes_final"
spe_celltype_path <- file.path(phenotype_dir, paste0(spe_celltype_name, ".rds"))
spe <- readRDS(spe_celltype_path)

metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_path <- file.path(metadata_dir, "ROI_Info_Aggregation.csv")
roi_df <- read_csv(roi_info_path)

# ROI list
roi_lst <- roi_df$ROI_ID

# ROI stage information
roi_vec <- c()
area_vec <- c()
stage_vec <- c()
for (ind in 1:nrow(roi_df)) {
    if (roi_df$ROI_Location[ind] == "Normal" | roi_df$ROI_Location[ind] == "DistantNormal") {
        roi_vec <- append(roi_vec, roi_df$ROI_ID[ind])
        area_vec <- append(area_vec, roi_df$Area[ind])
        stage_vec <- append(stage_vec, "Normal")        
    }
    else if (roi_df$ROI_Location[ind] == "Tumor") {
        roi_vec <- append(roi_vec, roi_df$ROI_ID[ind])
        area_vec <- append(area_vec, roi_df$Area[ind])
        stage_vec <- append(stage_vec, roi_df$ROI_Diag[ind])
    }
}

## obtain immune ratio
cell_id_lst <- rownames(colData(spe))
celtype_lst <- spe$celltype
immune_lst <- c("Macrophage", "Monocyte", "Dendritic-Cell", "Neutrophil", "MDSC", "NK-Cell")

immune_desnity_vec <- c()
for (ir in 1:length(roi_vec)) {
    cur_roi = roi_vec[ir]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- celtype_lst[cell_indices]
    immune_num <- 0
    for (immune_type in immune_lst) 
        immune_num <- immune_num + sum(roi_celltypes == immune_type)
    immune_desnity_vec <- append(immune_desnity_vec, immune_num * 1000000 / area_vec[ir])
}

# Construct data frame
stage_immune_density_df <- data.frame(ROI = roi_vec, Stage = stage_vec, Density = immune_desnity_vec)

MinMeanSEMMax <- function(x) {
    v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
}

stage_order <- c("Normal", "AAH", "AIS", "MIA", "ADC")
ggplot(stage_immune_density_df, aes(x = factor(Stage, level=stage_order), y=Density)) + 
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="black") + 
    geom_beeswarm(cex = 0.8, size=1.6)   

immune_density_plot_path <- file.path(data_root_dir, "NatureFigures", "Fig01", "Stage-Innate-Density.pdf")
ggsave(filename = immune_density_plot_path, device='pdf', width=8, height=8, dpi=300)
while (!is.null(dev.list()))
    dev.off()


immune_density_pfile <- file(file.path(data_root_dir, "NatureFigures", "Fig01", "Stage-Innate-Density-pval.txt"), "w")
# Welch Two Sample t-test p-value
normal_aah_ttest <- t.test(stage_immune_density_df$Density[stage_immune_density_df$Stage=="Normal"], stage_immune_density_df$Density[stage_immune_density_df$Stage=="AAH"])
writeLines(paste0("Normal-AAH:", normal_aah_ttest$p.value), immune_density_pfile)
normal_ais_ttest <- t.test(stage_immune_density_df$Density[stage_immune_density_df$Stage=="Normal"], stage_immune_density_df$Density[stage_immune_density_df$Stage=="AIS"])
writeLines(paste0("Normal-AIS:", normal_ais_ttest$p.value), immune_density_pfile)
normal_mia_ttest <- t.test(stage_immune_density_df$Density[stage_immune_density_df$Stage=="Normal"], stage_immune_density_df$Density[stage_immune_density_df$Stage=="MIA"])
writeLines(paste0("Normal-MIA:", normal_mia_ttest$p.value), immune_density_pfile)
normal_adc_ttest <- t.test(stage_immune_density_df$Density[stage_immune_density_df$Stage=="Normal"], stage_immune_density_df$Density[stage_immune_density_df$Stage=="ADC"])
writeLines(paste0("Normal-ADC:", normal_adc_ttest$p.value), immune_density_pfile)
aah_ais_ttest <- t.test(stage_immune_density_df$Density[stage_immune_density_df$Stage=="AAH"], stage_immune_density_df$Density[stage_immune_density_df$Stage=="AIS"])
writeLines(paste0("AAH-AIS:", aah_ais_ttest$p.value), immune_density_pfile)
aah_mia_ttest <- t.test(stage_immune_density_df$Density[stage_immune_density_df$Stage=="AAH"], stage_immune_density_df$Density[stage_immune_density_df$Stage=="MIA"])
writeLines(paste0("AAH-MIA:", aah_mia_ttest$p.value), immune_density_pfile)
aah_adc_ttest <- t.test(stage_immune_density_df$Density[stage_immune_density_df$Stage=="AAH"], stage_immune_density_df$Density[stage_immune_density_df$Stage=="ADC"])
writeLines(paste0("AAH-ADC:", aah_adc_ttest$p.value), immune_density_pfile)
ais_mia_ttest <- t.test(stage_immune_density_df$Density[stage_immune_density_df$Stage=="AIS"], stage_immune_density_df$Density[stage_immune_density_df$Stage=="MIA"])
writeLines(paste0("AIS-MIA:", ais_mia_ttest$p.value), immune_density_pfile)
ais_adc_ttest <- t.test(stage_immune_density_df$Density[stage_immune_density_df$Stage=="AIS"], stage_immune_density_df$Density[stage_immune_density_df$Stage=="ADC"])
writeLines(paste0("AIS-ADC:", ais_adc_ttest$p.value), immune_density_pfile)
mia_adc_ttest <- t.test(stage_immune_density_df$Density[stage_immune_density_df$Stage=="MIA"], stage_immune_density_df$Density[stage_immune_density_df$Stage=="ADC"])
writeLines(paste0("MIA-ADC:", mia_adc_ttest$p.value), immune_density_pfile)
close(immune_density_pfile)