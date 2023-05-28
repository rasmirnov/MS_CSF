install.packages("reticulate")
# install_github("ctlab/fgsea")
library(fgsea)
library(devtools)
library(reticulate)
library(Seurat)
library(anndata)
library(dplyr)
library(tidyverse)
library(viridis)
library(ggplot2)
library(ggrepel)
cluster <- 'myeloid_cells'

path <- "/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/romans/csf_atlas_new/"
setwd(paste0(path, cluster, '/whole/'))

#---------------------------------------1. DEA and HEATMAP---------------------------------------------------------
# upload raw h5ad obj
adata <- read_h5ad(paste0(cluster, "_raw.h5ad"))
# 1.create seurat obj with raw counts
seurat_data <- CreateSeuratObject(counts = t(as.matrix(adata$X)),      
                                  meta.data = adata$obs)

# # ## !!! optional: replace cDC2_2 to cDC2
seurat_data$cell_type <- sub("^(cDC2_2).*", "cDC2", seurat_data$cell_type)
# # ### !!! optional: filter out Neutrophils
seurat_data@meta.data <- seurat_data@meta.data %>%
  filter((cell_type == 'cDC2'))            # !(cell_type == 'Neutrophil')
# add column cell_type.organ here:
seurat_data@meta.data$cell_type.organ <- paste0(seurat_data@meta.data$cell_type, '_', 
                                                seurat_data@meta.data$organ)
# levels(pbulk) <- c("cDC2_PBMC", "cDC2_CSF")

# 2.create pseudobulk
pbulk <- AggregateExpression(
  object = seurat_data,
  assays = 'RNA',
  features = NULL,
  return.seurat = TRUE,
  group.by = c('cell_type.organ', 'sample'),
  add.ident = NULL,              
  slot = "counts",
  verbose = TRUE)

levels(pbulk)
Idents(pbulk)

# pbulk@meta.data %>% View()
pbulk$RNA@counts %>% max()
# round counts (for Deseq2 only):
pbulk$RNA@counts <- round(pbulk$RNA@counts, digits = 0)
pbulk$RNA@data %>% max()

# 3.normilize (for DE)
pbulk <- NormalizeData(pbulk,
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)

# 4.scale (for heatmap)
all.genes <- rownames(pbulk)
pbulk <- ScaleData(pbulk, features = all.genes)

# pbulk@assays$RNA@scale.data[1:10, 1:10]
pbulk@meta.data$cell_type.organ.sample <- rownames(pbulk@meta.data)
pbulk@meta.data$cell_type.organ <- sapply(pbulk@meta.data$cell_type.organ.sample,
                                          function(x) {str_extract(x, "[^_]*_[^_]*")})        # gsub("_[^_]+$", "\\1", x)

# 5.reassign Idents to "cell_type.organ" (for DE: comparation)
Idents(pbulk) <- pbulk[['cell_type.organ']]
levels(pbulk) <- c("cDC2_PBMC", "cDC2_CSF")
levels(pbulk)

# # 6.DE
# for MAST only (Bug): rename cell names
cnames <- paste0('test_', rownames(pbulk@meta.data))
pbulk <- RenameCells(pbulk, new.names = cnames)

# MAST = 1363 deg
# Deseq2 = 1676 deg
# how to get avg_log2FC with neg.values?
pbulk.markers <- FindAllMarkers(pbulk, 
                                test.use = "DESeq2",         # for DESeq2 slot is counts + rounded conts
                                slot = "slot",               # for MAST slot is data
                                only.pos = FALSE,
                                min.pct = 0.1,
                                logfc.threshold = 0.25)

pbulk.markers %>% View()

# 7.create 2 separate colms: cell_type & organ (for Phantasus/pheatmap)
pbulk.markers <- pbulk.markers %>% 
  mutate(organ = ifelse(grepl('PBMC', cluster), 'PBMC', 'CSF')) %>% 
  mutate(cell_type = sapply(cluster, function(x) {gsub("_.*", "", x)}))

pbulk.markers %>%
  group_by(organ, cell_type) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# 8.select top 25 from each group
top25 <- pbulk.markers %>%
  group_by(organ, cell_type) %>%
  top_n(n = 25, wt = avg_log2FC)

### check overal of markers between MAST and DESeq2 method
# # merge 
# top25_merged <-rbind(top25, top25_2)
# # find unique
# top25_unique <- top25_merged$gene %>% unique()
# # print them 
# cat(top25_unique, "\n")

levels(pbulk) <- c("cDC2_PBMC", "cDC2_CSF")
levels(pbulk)
pbulk@meta.data$cell_type.organ <- factor(pbulk@meta.data$cell_type.organ,
                                          levels = c("cDC2_PBMC", "cDC2_CSF"))
levels(pbulk@meta.data$cell_type.organ)


# 9.plot heatmap
DoHeatmap(pbulk, 
          features = top25$gene,
          group.by = "cell_type.organ",
          size = 4,       #8
          angle = 45,
          draw.lines = T) +
  scale_fill_viridis(option = 'inferno') +                    # option = 'inferno'
  guides(color = FALSE) +                     # hide cell_type legend
  theme(text = element_text(size = 17))

#---------------------------------------2. VOLCANO PLOT------------------------------------

# 1. Since they're symetrical & equal: filter out the half
pbulk.markers <- pbulk.markers %>% 
  filter(organ == 'CSF')


# # 2. centre DESeq2 data:
# pbulk.markers$avg_log2FC <- pbulk.markers$avg_log2FC - 0.6

# 3. Create new categorical column 
pbulk.markers <- pbulk.markers %>%
  mutate(gene_type = case_when(avg_log2FC >= 0 & p_val_adj <= 0.05 ~ "up in CSF",
                               avg_log2FC <= 0 & p_val_adj <= 0.05 ~ "up in PBMC",
                               TRUE ~ "p.adj > 0.05"))   

# # Statistics          
pbulk.markers %>%
  count(gene_type)

# 4. fix the order
pbulk.markers$gene_type <- factor(pbulk.markers$gene_type,
                                      levels = c("up in CSF", "up in PBMC", "p.adj > 0.05"))

# 5. define aeustetis parameters
cols <- c("up in CSF" = "#26b3ff", "up in PBMC" = "#ffad73", "p.adj > 0.05" = "grey") 
sizes <- c("up in CSF" = 2, "up in PBMC" = 2, "p.adj > 0.05" = 1) 
alphas <- c("up in CSF" = 0.8, "up in PBMC" = 0.8, "p.adj > 0.05" = 0.5)

#---------------------------------------optional-------------------------------------------------------
# 6. add labels of top N UP regulated genes in CSF and PBMC
top <- 10
top_genes <- bind_rows(
  pbulk.markers %>%
    filter(cluster == 'cDC2_CSF') %>%
    filter(gene_type == 'up in CSF') %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
    head(top),

  pbulk.markers %>%
    filter(cluster == 'cDC2_CSF') %>%
    filter(gene_type == 'up in PBMC') %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
    head(top)
)

# # # remove colms
# pbulk.markers <- pbulk.markers %>%
#   select(-c('top_genes', 'labels'))

# 7. Map those top N genes
# make a column with labeled rows (whether it's top gene or not)
pbulk.markers <- pbulk.markers %>%
  mutate(top_genes = case_when(rownames(pbulk.markers) %in% rownames(top_genes) ~ "top",   
                               TRUE ~ 'ns')) 

# create a col "labels" 
pbulk.markers$labels <- NA
pbulk.markers$labels[pbulk.markers$top_genes != "ns"] <- pbulk.markers$gene[pbulk.markers$top_genes != "ns"]

#----------------------------------------------------------------------------------------------

# 8. plot
pbulk.markers %>%
  ggplot(aes(x = avg_log2FC,       
             y = -log10(p_val_adj),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type,
             label = labels)) + 
  geom_point(shape = 21,    
             colour = "black") + 
  geom_text_repel(size = 4, box.padding = 1) +             # geom_label_repel()
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) +
  scale_x_continuous(breaks = c(seq(-6, 6, 2)),          # -10, 10, 2    
                     limits = c(-6, 6)) 


# 9. save in PNG/PDF
path_2 <- '/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/romans/figures/'
ggsave(paste0(path_2, "myeloid/volcano_DC2_CSF_deseq2_centered_2.png"),
       width = 15,    #12
       height = 8,    #6
       dpi = 300)

#--------------------------------------3. PATHWAY ANALYSIS (GSEA)----------------------------

de <- pbulk.markers 

stats <- de$p_val_adj
names(stats) <- de$gene

load("keggSymbolHuman.rdata")       # keggSymbolHuman.rdata
fgseaResults <- fgseaMultilevel(keggSymbolHuman, stats,
                                minSize = 1, maxSize = 500,
                                nPermSimple = 10000)
# head(fgseaResults, 3)

topPathwaysUp <- fgseaResults[ES > 0, ][head(order(pval), n=10), pathway]
# topPathwaysDown <- fgseaResults[ES < 0, ][head(order(pval), n=5), pathway]
topPathways <- topPathwaysUp        #c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(keggSymbolHuman[topPathways], stats, fgseaResults, gseaParam = 0.5)

fgseaResults_filtered <- fgseaResults %>% 
  filter(padj < 0.05)

ggplot(fgseaResults_filtered, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()


plotEnrichment(keggSymbolHuman[["Complement and coagulation cascades - Homo sapiens (human)"]],
               stats) + labs(title="Complement and coagulation cascades")


path_2 <- '/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/romans/figures/'
ggsave(paste0(path_2, "myeloid/fgsea_DC2.png"),
       width = 26,    #12
       height = 3,    #6
       dpi = 300)





















