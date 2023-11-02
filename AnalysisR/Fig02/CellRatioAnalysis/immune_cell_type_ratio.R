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
roi_vec <- roi_df$ROI_ID


## obtain lymphoid ratio
cell_id_lst <- rownames(colData(spe))
celtype_lst <- spe$celltype


immune_lst <- c("CD4-T-Cell", "CD8-T-Cell", "T-Reg-Cell", "B-Cell", "Macrophage", 
                "Monocyte", "Dendritic-Cell", "Neutrophil", "MDSC", "NK-Cell")
immune_type_lst <- c()
immune_ratio_lst <- c()

for (ir in 1:length(roi_vec)) {
    cur_roi = roi_vec[ir]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- celtype_lst[cell_indices]
    ttl_num <- 0.001
    for (cell_type in roi_celltypes) {
        if (cell_type %in% immune_lst)
            ttl_num <- ttl_num + 1
    }
    for (cell_type in immune_lst) {
        immune_type_lst <- append(immune_type_lst, cell_type)
        immune_ratio_lst <- append(immune_ratio_lst, sum(roi_celltypes == cell_type) * 1.0 / ttl_num)
    }
}

celltype_mean_ratio_lst <- list()
for (cell_type in immune_lst) {
    cell_indices <- which(immune_type_lst == cell_type)
    cell_ratio <- mean(immune_ratio_lst[cell_indices])
    celltype_mean_ratio_lst <- append(celltype_mean_ratio_lst, cell_ratio)
}
names(celltype_mean_ratio_lst) <- immune_lst
celltype_mean_ratio_lst <- celltype_mean_ratio_lst[order(unlist(celltype_mean_ratio_lst), decreasing=TRUE)]

# Construct data frame
cell_type_ratio_df <- data.frame(Type=immune_type_lst, Ratio=immune_ratio_lst)

for (cell_type in immune_lst) {
    print(cell_type)
    ct_df <- filter(cell_type_ratio_df, Type == cell_type)
    ct_ratios <- ct_df$Ratio
    ct_mean <- mean(ct_ratios)
    ct_sd <- sd(ct_ratios)
    print(paste0("Mean: ", ct_mean, "SD: ", ct_sd))
}

MinMeanSEMMax <- function(x) {
    v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
}

ggplot(cell_type_ratio_df, aes(x = factor(Type, level=names(celltype_mean_ratio_lst)), y=Ratio)) + 
    theme(axis.text.x=element_text(angle=90, hjust=1)) + 
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="black") + 
    geom_beeswarm(cex = 2.5, corral = "random", corral.width = 0.4)