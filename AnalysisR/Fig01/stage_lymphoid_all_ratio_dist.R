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
roi_info_name <- "ROI_Info"
roi_info_path <- file.path(metadata_dir, paste0(roi_info_name, ".xlsx"))
roi_df <- read_excel(roi_info_path)

# ROI list
roi_lst <- roi_df$ROI_ID

# ROI stage information
roi_vec <- c()
stage_vec <- c()
for (ind in 1:nrow(roi_df)) {
    if (roi_df$ROI_Location[ind] == "Normal" | roi_df$ROI_Location[ind] == "DistantNormal") {
        roi_vec <- append(roi_vec, roi_df$ROI_ID[ind])
        stage_vec <- append(stage_vec, "Normal")        
    }
    else if (roi_df$ROI_Location[ind] == "Tumor") {
        roi_vec <- append(roi_vec, roi_df$ROI_ID[ind])
        stage_vec <- append(stage_vec, roi_df$ROI_Diag[ind])
    }
}

## obtain lymphoid ratio
cell_id_lst <- rownames(colData(spe))
celtype_lst <- spe$celltype
lymphoid_lst <- c("CD4-T-Cell", "CD8-T-Cell", "T-Reg-Cell", "B-Cell", "NK-Cell")
lymphoid_ratio_vec <- c()
for (ir in 1:length(roi_vec)) {
    cur_roi = roi_vec[ir]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- celtype_lst[cell_indices]
    ttl_num = length(roi_celltypes)
    lymphoid_num <- 0
    for (lymphoid_type in lymphoid_lst) 
        lymphoid_num <- lymphoid_num + sum(roi_celltypes == lymphoid_type)
    lymphoid_ratio_vec <- append(lymphoid_ratio_vec, lymphoid_num * 1.0 / ttl_num)
}

# Construct data frame
stage_lymphoid_ratio_df <- data.frame(ROI = roi_vec, Stage = stage_vec, Ratio = lymphoid_ratio_vec)

MinMeanSEMMax <- function(x) {
    v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
}

stage_order <- c("Normal", "AAH", "AIS", "MIA", "ADC")
ggplot(stage_lymphoid_ratio_df, aes(x = factor(Stage, level=stage_order), y=Ratio)) + 
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="black") + 
    geom_beeswarm(cex = 0.8, size=1.6)   

immune_ratio_plot_path <- file.path(data_root_dir, "NatureFigures", "Fig01", "Stage-Lymphoid-All-Ratio.pdf")
ggsave(filename = immune_ratio_plot_path, device='pdf', width=8, height=8, dpi=300)
while (!is.null(dev.list()))
    dev.off()


immune_ratio_pfile <- file(file.path(data_root_dir, "NatureFigures", "Fig01", "Stage-Lymphoid-All-Ratio-pval.txt"), "w")
# Welch Two Sample t-test p-value
normal_aah_ttest <- t.test(stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="Normal"], stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="AAH"])
writeLines(paste0("Normal-AAH:", normal_aah_ttest$p.value), immune_ratio_pfile)
normal_ais_ttest <- t.test(stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="Normal"], stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="AIS"])
writeLines(paste0("Normal-AIS:", normal_ais_ttest$p.value), immune_ratio_pfile)
normal_mia_ttest <- t.test(stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="Normal"], stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="MIA"])
writeLines(paste0("Normal-MIA:", normal_mia_ttest$p.value), immune_ratio_pfile)
normal_adc_ttest <- t.test(stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="Normal"], stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="ADC"])
writeLines(paste0("Normal-ADC:", normal_adc_ttest$p.value), immune_ratio_pfile)
aah_ais_ttest <- t.test(stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="AAH"], stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="AIS"])
writeLines(paste0("AAH-AIS:", aah_ais_ttest$p.value), immune_ratio_pfile)
aah_mia_ttest <- t.test(stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="AAH"], stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="MIA"])
writeLines(paste0("AAH-MIA:", aah_mia_ttest$p.value), immune_ratio_pfile)
aah_adc_ttest <- t.test(stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="AAH"], stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="ADC"])
writeLines(paste0("AAH-ADC:", aah_adc_ttest$p.value), immune_ratio_pfile)
ais_mia_ttest <- t.test(stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="AIS"], stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="MIA"])
writeLines(paste0("AIS-MIA:", ais_mia_ttest$p.value), immune_ratio_pfile)
ais_adc_ttest <- t.test(stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="AIS"], stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="ADC"])
writeLines(paste0("AIS-ADC:", ais_adc_ttest$p.value), immune_ratio_pfile)
mia_adc_ttest <- t.test(stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="MIA"], stage_lymphoid_ratio_df$Ratio[stage_lymphoid_ratio_df$Stage=="ADC"])
writeLines(paste0("MIA-ADC:", mia_adc_ttest$p.value), immune_ratio_pfile)
close(immune_ratio_pfile)

