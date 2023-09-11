library(imcRtools)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

## set directory
data_root_dir <- "E:/LungIMCData/HumanWholeIMC"
phenotype_dir <- file.path(data_root_dir, "CellPhenotyping")
celltype_delaunay_dir <- file.path(phenotype_dir, "DelaunayInteraction")
# load the spe
spe_cn_path <- file.path(celltype_delaunay_dir, "Delaunay50-CN8.rds")
spe <- readRDS(spe_cn_path)

metadata_dir <- file.path(data_root_dir, "Metadata")
roi_info_name <- "ROI_Info"
roi_info_path <- file.path(metadata_dir, paste0(roi_info_name, ".xlsx"))
roi_df <- read_excel(roi_info_path)
# ROI list
roi_vec <- roi_df$ROI_ID
## obtain cell type list
cell_id_lst <- rownames(colData(spe))
cn_type_lst <- spe$cn_celltypes


all_cell_lst <- c(1, 2, 3, 4, 5, 6, 7, 8)
all_type_lst <- c()
all_ratio_lst <- c()

for (ir in 1:length(roi_vec)) {
    cur_roi = roi_vec[ir]
    cell_indices <- which(startsWith(cell_id_lst, cur_roi))
    roi_celltypes <- cn_type_lst[cell_indices]
    ttl_num = length(roi_celltypes)
    for (cell_type in all_cell_lst) {
        all_type_lst <- append(all_type_lst, cell_type)
        all_ratio_lst <- append(all_ratio_lst, sum(roi_celltypes == cell_type) * 1.0 / ttl_num)
    }
}

# update type name
all_type_lst[all_type_lst == 1] <- "Undefined-CN1"
all_type_lst[all_type_lst == 2] <- "Epithelial1-CN2"
all_type_lst[all_type_lst == 3] <- "Proliferating-CN3"
all_type_lst[all_type_lst == 4] <- "Epithelial2-CN4"
all_type_lst[all_type_lst == 5] <- "Endothelial-CN5"
all_type_lst[all_type_lst == 6] <- "Fibroblast-CN6"
all_type_lst[all_type_lst == 7] <- "Macrophage-CN7"
all_type_lst[all_type_lst == 8] <- "Pan-Immune-CN8"

all_cell_lst <- c("Undefined-CN1", "Epithelial1-CN2", "Proliferating-CN3", "Epithelial2-CN4",
                  "Endothelial-CN5", "Fibroblast-CN6", "Macrophage-CN7", "Pan-Immune-CN8")

print("Proportion mean & std:")
celltype_mean_ratio_lst <- list()
for (cell_type in all_cell_lst) {
    cell_indices <- which(all_type_lst == cell_type)
    cell_ratio_mean <- mean(all_ratio_lst[cell_indices])
    cell_ratio_std <- sd(all_ratio_lst[cell_indices])
    print(paste(cell_type, ": mean", cell_ratio_mean, "std", cell_ratio_std))
    celltype_mean_ratio_lst <- append(celltype_mean_ratio_lst, cell_ratio_mean)
}
names(celltype_mean_ratio_lst) <- all_cell_lst
celltype_mean_ratio_lst <- celltype_mean_ratio_lst[order(unlist(celltype_mean_ratio_lst), decreasing=TRUE)]


# Construct data frame
cell_type_ratio_df <- data.frame(Type=all_type_lst, Ratio=all_ratio_lst)
MinMeanSEMMax <- function(x) {
    v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
}

ggplot(cell_type_ratio_df, aes(x = factor(Type, level=names(celltype_mean_ratio_lst)), y=Ratio)) + 
    theme(axis.text.x=element_text(angle=90, hjust=1)) + 
    stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", colour="black") + 
    geom_beeswarm(cex = 2.5, corral = "random", corral.width = 0.4)
