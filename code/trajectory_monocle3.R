# install.packages('textshaping', type = 'source')
# devtools::install_version('speedglm', '0.3-4',
#                           repos = 'https://packagemanager.rstudio.com/cran/2023-03-31')
# install.packages("/home/roma/Downloads/satijalab-seurat-wrappers-d28512f.tar.gz", 
#                  repos = NULL, type = "source")
setwd('~/Downloads/MD_proj/neuro_proj/data/data/micro_mono_csf/')
library(BiocGenerics)
library(DelayedArray)
library(DelayedMatrixStats)
library(limma)
library(S4Vectors)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(HDF5Array)
library(terra)     
library(ggrastr)
library(devtools)
library(speedglm)
library(monocle3)
#-----------------
library(ggplot2)
library(Seurat)
library(dplyr)
library(reticulate)
library(anndata)
library(harmony)
library(R.utils)
library(SeuratWrappers)
#----------------
# 1.export micro_csf_mono & embeddings from server:
# clusters <-  c('micro_mono_csf')
# for (cluster in clusters) {
  # change dir to merged
  path <- "~/Downloads/MD_proj/neuro_proj/data/data/micro_mono_csf/"
  setwd(path)
  
  adata <- read_h5ad("~/Downloads/MD_proj/neuro_proj/data/data/micro_mono_csf/micro_mono_csf.h5ad")
  
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
  
  saveRDS(seurat_data, paste0(cluster, ".rds"))
}

seurat_data <- readRDS('micro_mono_csf.rds')

# plot UMAP based on harmony embeddings:
seurat_data <- RunUMAP(seurat_data, n.neighbors = 10L, dims = 1:30, reduction = 'pca_harmony')

DimPlot(seurat_data, pt.size = 0.8, reduction = "umap", group.by = 'cell_type')

seurat_data@active.assay <- 'RNA'

# 4.convert to monocle3 obj
cds <- SeuratWrappers::as.cell_data_set(seurat_data)    # meta stored here: cds@colData$cell_type

# cluster:
cds <- cluster_cells(cds)

plot_cells(cds, show_trajectory_graph = FALSE, 
           color_cells_by = 'cell_type',
           group_label_size = 4)

# create a graph:
cds <- learn_graph(cds, use_partition = F)

# set the root cell type:
# cds <- order_cells(cds)            # choose by interface
cds <- order_cells(cds, root_cells = colnames(cds[ ,clusters(cds) == 3]))   # assign CD14 Mono as a root

# plot pseudotime trajectory:
plot_cells(cds, color_cells_by = 'pseudotime', 
           label_branch_points = F, 
           label_leaves = T, 
           label_roots = F,
           cell_size = 0.8,
           trajectory_graph_segment_size = 1,
           graph_label_size = 4)

# calculate pseudotime statistics:
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds) 
df.psedo <- as.data.frame(colData(cds))

# plot pseudotime statistics per cell type:
ggplot(df.psedo, aes(reorder(cell_type, monocle3_pseudotime, median), monocle3_pseudotime, fill = cell_type)) + 
  geom_boxplot() + 
  geom_point(aes(x = cell_type), shape = 21, position = position_dodge(width = 1)) + 
  xlab('Cell type') +
  ylab('Pseudotime')

# save
ggsave(paste0("~/Downloads/MD_proj/neuro_proj/results/figures/trajectory/micro_ptime_leaves.png"),
       width = 10,    #10
       height = 7,    #7
       dpi = 300)
