library(ggplot2)
library(Seurat)
library(dplyr)
library(reticulate)
library(anndata)
library(tidyverse)
library(devtools)
library(dyno)

library(hdf5r)
library(hexbin)

# 1.upload RDS obj:
cluster <-  c('micro_bam_csf_cd14')
path <- "~/Downloads/MD_proj/neuro_proj/data/data/micro_mono_csf/"
setwd(path)

seurat_data <- readRDS(paste0(cluster, '_raw.rds'))

seurat_data@reductions$umap@key <- 'UMAP_'
colnames(x = seurat_data[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)

seurat_data@meta.data$cell_type <- as.character(seurat_data@meta.data$cell_type)

seurat_data@active.assay <- 'RNA'

### check the position of chosen cell name on the UMAP----------------------------

# cell_names <- 'CGCGTGCGCTTA_2-0-0-0'           #ATGCTAGTTATC_4-0-0-1
# 
# DimPlot(seurat_data,
#         cells.highlight = cell_names,
#         pt.size = 2,
#         cols.highlight = "red")
#-------------------------------------------------------------------------------------------
# plot UMAP based on harmony embeddings:
# DimPlot(seurat_data,
#         pt.size = 0.8,
#         reduction = "umap",
#         group.by = 'cell_type')

# 2. Normilize RAW data
seurat_data <- NormalizeData(seurat_data,
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)

# 3. Creat wrap_data obj for dyno
Idents(seurat_data) <- seurat_data@meta.data$cell_type

dataset <- wrap_expression(
  counts = t(as.matrix(seurat_data@assays$RNA@counts)),
  expression = t(as.matrix(seurat_data@assays$RNA@data))
)

dataset <- add_prior_information(dataset, 
                                 start_id = 'CGCGTGCGCTTA_2-0-0-0',         # CSF_Mono
                                 # sample(colnames(subset(seurat_data,
                                 #                                   idents = 'CD14 Mono CSF')), 1),
                                 end_id = 'CCGTTACAATCC_5-0-0-1',           # SPP1+
                                 # sample(colnames(subset(seurat_data,
                                 #                                 idents = 'SPP1+ MG')), 1),
                                 start_n = 1,
                                 end_n = 1,
                                 groups_id=cbind(colnames(seurat_data), seurat_data[['cell_type']][[1]]) %>%
                                   as.data.frame() %>%
                                   magrittr::set_colnames(c('cell_id', 'group_id')))
dataset <- add_dimred(
  dataset,
  dimred = as.matrix(seurat_data@reductions$umap@cell.embeddings))

dataset <- add_grouping(
  dataset,
  seurat_data[['cell_type']][[1]]
)

# ----------------------optional-----------------------------------------------
# guidelines <- guidelines_shiny(dataset)
# methods_selected <- guidelines$methods_selected[1]    ### click on the bottom "use/continue" 
# dynwrap::test_docker_installation(detailed = TRUE)
#--------------------------------------------------------------------------------

# 4. Run trajectory inference analysis
model <- infer_trajectory(dataset, 
                          'slingshot',                    # try monocle_ddrtree ?
                          # first(methods_selected),     
                          give_priors = c('start_id', 'end_id',
                                          'start_n', 'end_n'))

# ----------------------optional-----------------------------------------------
# add_root
add_root(
  model,
  root_cell_id = model$root_cell_id,                  
  # root_milestone_id = model$root_milestone_id,
  flip_edges = TRUE
)
# -----------------------------------------------------------------------------

# 5.1 plot TI
plot_dimred(model,
            dimred = as.matrix(seurat_data@reductions$umap@cell.embeddings),
            expression_source = dataset$expression,     # normilized
            grouping = seurat_data$cell_type,
            color_density = "grouping",
            alpha_cells = 0.7,
            size_cells = 1,
            label_milestones = F,
            hex_cells = F)

# 5.2 plot Pseudotime
patchwork::wrap_plots(
  plot_dimred(model, 
              dimred = as.matrix(seurat_data@reductions$umap@cell.embeddings),
              "pseudotime", 
              pseudotime = calculate_pseudotime(model),
              alpha_cells = 0.7,
              size_cells = 1,
              label_milestones = F,
              hex_cells = F) + 
    ggtitle("Pseudotime")
)

# 6. save
ggsave("~/Downloads/MD_proj/neuro_proj/results/figures/by_object/micro/trajectory/micro_bam_cd14_csf/pdf/slingshot.pdf",
       width = 9,    #10
       height = 6,    #7
       dpi = 300)

