# load Panel metadata
panel_file_path <- file.path(data_root, "BatchCorrection", "Metadata", "PanelsIMC.csv")
markers <- read.csv(panel_file_path) %>% pull(Marker)