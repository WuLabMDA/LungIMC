fea_root_dir <- file.path(data_root, "BatchCorrection", "NoControl")
load(file.path(fea_root_dir, "CorrectedFeas.RData"))
community_name = "CorrectedCommunitiesPG.RData"
cell_feas <- select(corrected, CK, aSMA, CD31, CD45)
cell_feas <- as.matrix(cell_feas)
message(paste("Cell Number: ", NROW(cell_feas), " Feature Number: ", NCOL(cell_feas)))

# Setup FastPG parameters
k <- 30
num_threads <- 32

# Run FastPG
message(paste("FastPG Start @ ", format(Sys.time(), "%a %b %d %X %Y")))
clusters <- FastPG::fastCluster(cell_feas, k, num_threads)
message(paste("FastPG Finish @ ", format(Sys.time(), "%a %b %d %X %Y")))

# Analyze results
transform_communities <- clusters$communities
unique_ids <- unique(transform_communities)

if(all(unique_ids >= 0)) {
    message("No singleton clusters.")
} else {
    message("Singleton clusters exists!")
}

message(paste("There are ", length(unique_ids), " communties detected."))
community_path <- file.path(fea_root_dir, community_name)
save(cell_feas, transform_communities, file = community_path)
