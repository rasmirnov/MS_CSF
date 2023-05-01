clusters = c('myeloid')
library(Seurat)
library(anndata)
library(dplyr)

for (cluster in clusters) {
# change dir to merged
  path <- "/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/romans/csf_atlas_new/"
  setwd(paste0(path, cluster, '_cells/whole'))
  # setwd(paste0(path, 'csf_atlas'))

  adata <- read_h5ad(paste0(cluster, "_raw.h5ad"))
  
  pca <- as.matrix(read.csv("embeddings/pca.csv", row.names = 1))

  pca_harmony <- as.matrix(read.csv("embeddings/pca_harmony.csv", row.names = 1))

  umap <- as.matrix(read.csv("embeddings/umap.csv", row.names = 1))

  seurat_data <- CreateSeuratObject(counts = t(as.matrix(adata$X)),      # as.matrix
                                    meta.data = adata$obs)
  
  seurat_data[["pca"]] <- CreateDimReducObject(embeddings = pca,
                                                         key = "pca_", 
                                                         assay = DefaultAssay(seurat_data))
  
  seurat_data[["pca_harmony"]] <- CreateDimReducObject(embeddings = pca_harmony,
                                                         key = "pca_harmony_", 
                                                         assay = DefaultAssay(seurat_data))
  
  seurat_data[["umap"]] <- CreateDimReducObject(embeddings = umap,
                                                         key = "umap_", 
                                                         assay = DefaultAssay(seurat_data))
  
  colnames(seurat_data@meta.data) <- gsub('leiden', 'leiden_res', colnames(seurat_data@meta.data))
  
  seurat_data <- NormalizeData(seurat_data,
                normalization.method = "LogNormalize",
                scale.factor = 10000)
  # delete res 1.4-2.0 in metadata:
  seurat_data$leiden_res_0.2 <- NULL
  seurat_data$leiden_res_0.3 <- NULL
  seurat_data$leiden_res_0.5 <- NULL
  seurat_data$leiden_res_0.6 <- NULL
  # seurat_data$leiden_res_1.0 <- NULL
  # seurat_data$leiden_res_1.8 <- NULL
  # seurat_data$leiden_res_2.0 <- NULL
  
  saveRDS(seurat_data, paste0(cluster, ".rds"))
}
# check resolutions & # cells