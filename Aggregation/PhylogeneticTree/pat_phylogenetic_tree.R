root_path = "E:/LungIMCData/HumanWholeIMC/Aggregation/Phylo"
source("tree_utils.R")


# pid <- "2405"
# pid <- "H16-0223"
pid <- "H18-0271"
fea_path <- file.path(root_path, paste(pid, "roi_feas.csv", sep="_"))
p_feas <- read.table(fea_path, header=TRUE, sep=",", quote="", stringsAsFactors=FALSE)
col_num <- dim(p_feas)[2]

# reconstruct feature names
roi_names = c()
for (ind in 1:length(p_feas$ROI_ID)) {
    p_roi_len = nchar(p_feas$ROI_ID[ind])
    p_roi_name <- substr(p_feas$ROI_ID[ind], nchar(pid)+2, p_roi_len)
    roi_names <- append(roi_names, paste(p_roi_name, p_feas$ROI_Stage[ind], sep="-"))
}
rownames(p_feas) = roi_names
p_feas = p_feas[, 3:col_num]

tree_path <- file.path(root_path, paste(pid, "phylogenetic_tree.pdf", sep="_"))
pdf(tree_path, height=15, width=18)
pretty_tree(p_feas, num_clusters=3)
dev.off()