# load data
rdata_dir <- file.path(data_root, "RData")
no_control_meta_path <- file.path(rdata_dir, "StudyMeta.RData")
load(no_control_meta_path)
raw_fea_path <- file.path(rdata_dir, "TransformFeas.csv")
raw_feas <- read.csv(raw_fea_path, header=TRUE)
cell_feas <- select(raw_feas, all_of(markers))
communities <- raw_feas$label
message(paste("There are ", length(unique(communities)), " communties detected."))
fea_community_name = "RawFeaCommunities.RData"
fea_community_path <- file.path(rdata_dir, fea_community_name)
save(cell_feas, communities, file = fea_community_path)