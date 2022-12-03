library(imcRtools)
library(cytomapper)
library(stringr)
library(openxlsx)


data_root_dir <- "E:/LungIMCData/HumanSampling35-0/LungROIProcessing/Steinbock"


# load antibody information
spe_file_path <- file.path(data_root_dir, "steinbock_spe.rds")
spe <- readRDS(spe_file_path)

# load 
expr_vals <- assay(spe, "exprs")
antibody_lst <- rownames(expr_vals)
cell_lst <- colnames(expr_vals)

CD8a_expr_vals <- expr_vals[which(antibody_lst == "CD8a"),1:length(cell_lst)] 
CD8Tcells <- which(CD8a_expr_vals > 0.04) # 0.06
CD8Tcell_ids <- names(CD8Tcells)
CD8Tcell_indices <- as.numeric(CD8Tcells)

CD8Tcells_df <- data.frame(unlist(CD8Tcell_ids), unlist(CD8Tcell_indices))
names(CD8Tcells_df) <- c("IDS", "Index")
CD8T_file_path <- file.path(data_root_dir, "CD8Tcells.csv")
write.csv(CD8Tcells_df, CD8T_file_path, row.names = FALSE)


CD19_expr_vals <- expr_vals[which(antibody_lst == "CD19"),1:length(cell_lst)]
Bcells <- which(CD19_expr_vals > 0.02) # 0.02
Bcell_ids <- names(Bcells)
Bcell_indices <- as.numeric(Bcells)

common_cells <- intersect(CD8Tcell_indices, Bcell_indices)



