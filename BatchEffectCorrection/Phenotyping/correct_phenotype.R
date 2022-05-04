# load data
rdata_dir <- file.path(data_root, "RData")
no_control_meta_path <- file.path(rdata_dir, "StudyMeta.RData")
load(no_control_meta_path)
correct_fea_path <- file.path(rdata_dir, "SelfCorrectFeas.csv")
correct_feas <- read.csv(correct_fea_path, header=TRUE)
cell_feas <- select(correct_feas, all_of(markers))
communities <- correct_feas$label
message(paste("There are ", length(unique(communities)), " communties detected."))
fea_community_name = "CorrectFeaCommunities.RData"
fea_community_path <- file.path(rdata_dir, fea_community_name)
save(cell_feas, communities, file = fea_community_path)