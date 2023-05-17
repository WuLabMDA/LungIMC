root_path = "E:/LungIMCData/HumanWholeIMC/Aggregation/Phylo"
source("all_tree_utils.R")


fea_path <- file.path(root_path, "all_roi_feas.csv")
roi_feas <- read.table(fea_path, header=TRUE, sep=",", quote="", stringsAsFactors=FALSE)
col_num <- dim(roi_feas)[2]
roi_stages <- roi_feas[, 2]
roi_clusters <- rep(0, length(roi_stages))
for (m in 1:length(roi_stages)) {
    if (roi_stages[m] == "Normal")
        roi_clusters[m] <- 1
    else if (roi_stages[m] == "AAH")
        roi_clusters[m] <- 2
    else if (roi_stages[m] == "AIS")
        roi_clusters[m] <- 3    
    else if (roi_stages[m] == "MIA")
        roi_clusters[m] <- 4
    else if (roi_stages[m] == "ADC")
        roi_clusters[m] <- 5
    else
        print(roi_stages[m])
}
roi_feas <- roi_feas[, 3:col_num]

tree_path <- file.path(root_path, "all_phylogenetic_tree.pdf")
pdf(tree_path, height=15, width=18)


pretty_tree(roi_feas, roi_clusters)
dev.off()