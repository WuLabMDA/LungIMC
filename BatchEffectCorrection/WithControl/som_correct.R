rdata_dir <- file.path(data_root, "RData")
if (!dir.exists(rdata_dir))
    dir.create(rdata_dir, recursive = TRUE)

community_name = "ControlCorrectedCommunitiesSOM.RData"
cell_feas <- select(corrected, CK, aSMA, CD31, CD45)
cell_feas <- as.matrix(cell_feas)
message(paste("Cell Number: ", NROW(cell_feas), " Feature Number: ", NCOL(cell_feas)))
som_labels <- kohonen::som(cell_feas, grid = kohonen::somgrid(xdim = 8, ydim = 8), 
                           rlen = 10, dist.fcts = "euclidean")
correct_communities <- som_labels$unit.classif
message(paste("There are ", length(unique(correct_communities)), " communties detected."))
community_path <- file.path(rdata_dir, community_name)
save(cell_feas, correct_communities, file = community_path)
