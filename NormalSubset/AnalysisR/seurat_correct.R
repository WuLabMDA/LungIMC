library(patchwork)
library(Seurat)
library(SeuratObject)
library(scater)
library(cowplot)
library(dittoSeq)
library(viridis)
set.seed(1234)

normal_root_dir <- "E:/LungIMCData/LungROIProcessing/SteinbockNormal"
raw_spe_path = file.path(normal_root_dir, "steinbock_spe.rds")
spe <- readRDS(raw_spe_path)

# Convert to Seurat object
seurat_obj <- as.Seurat(spe, counts = "counts", data = "exprs")
seurat_obj <- AddMetaData(seurat_obj, as.data.frame(colData(spe)))
seurat.list <- SplitObject(seurat_obj, split.by = "stain_id")

# PCA -  finding anchors
features <- rownames(spe)[rowData(spe)$use_channel]
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    return(x)
})


print(paste("---Start find integration anchors---"))
anchors <- FindIntegrationAnchors(object.list = seurat.list,
                                  anchor.features = features,
                                  k.anchor = 8,
                                  k.filter = 64,
                                  k.score = 16,
                                  n.trees = 20)

print(paste("---Start integration---"))
combined <- IntegrateData(anchorset = anchors)
print(paste("Finish at", format(Sys.time(), "%a %b %d %X %Y")))

# Select the integrated assay and perform PCA dimensionality reduction
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 24, verbose = FALSE)
reducedDim(spe, "seurat") <- combined@reductions$pca@cell.embeddings

# Compute the UMAP embeddings based on Seurat integrated results
spe <- runUMAP(spe, dimred = "seurat", name = "UMAP_seurat")
seurat_spe_path = file.path(normal_root_dir, "seurat_spe.rds")
saveRDS(spe, seurat_spe_path)