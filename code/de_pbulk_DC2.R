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
# seurat_data$cell_type <- sub("^(cDC2_2).*", "cDC2", seurat_data$cell_type)
# # ### !!! optional: filter out Neutrophils
seurat_data@meta.data <- seurat_data@meta.data %>%
  filter((cell_type == c('CD32B+ cDC2 like', 'CD36+ cDC2')))            # !(cell_type == 'Neutrophil')
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
  group.by = c('cell_type', 'sample'),
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
pbulk@meta.data$cell_type.sample <- rownames(pbulk@meta.data)
pbulk@meta.data$cell_type <- sapply(pbulk@meta.data$cell_type.sample,
                                          function(x) {str_extract(x, "[^_]*")})        # gsub("_[^_]+$", "\\1", x)

# 5.reassign Idents to "cell_type.organ" (for DE: comparation)
Idents(pbulk) <- pbulk[['cell_type']]
levels(pbulk) <- c("CD36+ cDC2", "CD32B+ cDC2 like")
levels(pbulk)

# # 6.DE
# for MAST only (Bug): rename cell names
cnames <- paste0('test_', rownames(pbulk@meta.data))
pbulk <- RenameCells(pbulk, new.names = cnames)

# MAST = 1363 deg
# Deseq2 = 1676 deg
# how to get avg_log2FC with neg.values?
pbulk.markers <- FindAllMarkers(pbulk, 
                                test.use = "MAST",         # for DESeq2 slot is counts + rounded conts
                                slot = "data",               # for MAST slot is data
                                only.pos = FALSE,
                                min.pct = 0.1,
                                logfc.threshold = 0.25)

# pbulk.markers.deseq2 <- FindAllMarkers(pbulk, 
#                                 test.use = "DESeq2",         # for DESeq2 slot is counts + rounded conts
#                                 slot = "counts",               # for MAST slot is data
#                                 only.pos = FALSE,
#                                 min.pct = 0.1,
#                                 logfc.threshold = 0.25)

# pbulk.markers %>% View()
pbulk.markers %>% nrow()

# 7.create 2 separate colms: cell_type & organ (for Phantasus/pheatmap)
pbulk.markers <- pbulk.markers %>% 
  # mutate(organ = ifelse(grepl('PBMC', cluster), 'PBMC', 'CSF')) %>% 
  # mutate(cell_type = sapply(cluster, function(x) {gsub("_.*", "", x)})) %>% 
  rename(cell_type = cluster) 
#   # filter(p_val_adj < 0.05)
# pbulk.markers.deseq2 <- pbulk.markers.deseq2 %>% 
#   # mutate(organ = ifelse(grepl('PBMC', cluster), 'PBMC', 'CSF')) %>% 
#   # mutate(cell_type = sapply(cluster, function(x) {gsub("_.*", "", x)})) %>% 
#   rename(cell_type = cluster)

pbulk.markers %>%
  group_by(cell_type) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# 8.select top 25 from each group
top25 <- pbulk.markers %>%
  filter(p_val_adj < 0.05) %>% 
  group_by(cell_type) %>%
  top_n(n = 25, wt = avg_log2FC)

### check overal of markers between MAST and DESeq2 method
# # merge 
# top25_merged <-rbind(top25, top25_2)
# # find unique
# top25_unique <- top25_merged$gene %>% unique()
# # print them 
# cat(top25$gene, "\n")

levels(pbulk) <- c("CD36+ cDC2", "CD32B+ cDC2 like")
levels(pbulk)
pbulk@meta.data$cell_type <- factor(pbulk@meta.data$cell_type,
                                          levels = c("CD36+ cDC2", "CD32B+ cDC2 like"))
levels(pbulk@meta.data$cell_type)


# 9.plot heatmap
DoHeatmap(pbulk, 
          features = top25$gene,
          group.by = "cell_type",
          size = 8,       #8
          angle = 45,
          draw.lines = T) +
  scale_fill_viridis(option = 'inferno') +                    # option = 'inferno'
  guides(color = FALSE) +                     # hide cell_type legend
  theme(text = element_text(size = 10))

# 10. save in PNG/PDF
path_2 <- '/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/romans/figures/'
ggsave(paste0(path_2, "myeloid/whole/heatmap_DC2_MAST.png"),
       width = 20,    #12
       height = 12,    #6
       dpi = 300)
#---------------------------------------2. VOLCANO PLOT------------------------------------

# 1. Since they're symetrical & equal: filter out the half
pbulk.markers.new <- pbulk.markers %>%
  filter(cell_type == 'CD32B+ cDC2 like') 
  # filter(p_val_adj < 0.05)

# 3. Create new categorical column 
pbulk.markers.new <- pbulk.markers.new %>%
  mutate(gene_type = case_when(gene %in% c('APOC1', 'AM2',
                                           'RGS1', 'CST7', 'HPGDS', 'BHLHE41', 'A2M',
                                           'AXL', 'NR4A3', 'FRMD4A', 'FMNL3', 'BIN1') ~ "up in CD32B+ cDC2 like\n and coexpressed in Microglia",
                               avg_log2FC >= 0 & p_val_adj <= 0.05 ~ "up in CD32B+ cDC2 like",
                               avg_log2FC <= 0 & p_val_adj <= 0.05 ~ "up in CD36+ cDC2",
                               TRUE ~ "p.adj > 0.05"))   

# # Statistics          
pbulk.markers.new %>%
  count(gene_type)

# 4. fix the order
pbulk.markers.new$gene_type <- factor(pbulk.markers.new$gene_type,
                                      levels = c("up in CD36+ cDC2", "up in CD32B+ cDC2 like",
                                                 "up in CD32B+ cDC2 like\n and coexpressed in Microglia",
                                                 "p.adj > 0.05"))

# 5. define aeustetis parameters (#FFDB6D #00AFBB  ) 
cols <- c("up in CD32B+ cDC2 like" = "#26b3ff", "up in CD32B+ cDC2 like\n and coexpressed in Microglia" = '#ff7387',
          "up in CD36+ cDC2" = "#ffad73", "p.adj > 0.05" = "grey")
sizes <- c("up in CD32B+ cDC2 like" = 2, "up in CD32B+ cDC2 like\n and coexpressed in Microglia" = 2, 
           "up in CD36+ cDC2" = 2, "p.adj > 0.05" = 1) 
alphas <- c("up in CD32B+ cDC2 like" = 0.8, "up in CD32B+ cDC2 like\n and coexpressed in Microglia" = 0.8,
            "up in CD36+ cDC2" = 0.8, "p.adj > 0.05" = 0.5)

#---------------------------------------optional-------------------------------------------------------
# 6. add labels of top N UP regulated genes in CSF and PBMC
top <- 15
top_genes <- bind_rows(
  pbulk.markers.new %>%
    filter(gene_type == 'up in CD36+ cDC2') %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
    head(top), 

  pbulk.markers.new %>%
    filter(gene_type == 'up in CD32B+ cDC2 like') %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
    head(top),

  pbulk.markers.new %>%
    filter(gene_type == 'up in CD32B+ cDC2 like\n and coexpressed in Microglia')
)

# filter out black duplicates:
top_genes <-  top_genes[!duplicated(top_genes$gene),]

# # # remove colms
# pbulk.markers.new <- pbulk.markers.new %>%
#   select(-c('top_genes', 'labels'))

# 7. Map those top N genes
# make a column with labeled rows (whether it's top gene or not)
pbulk.markers.new <- pbulk.markers.new %>%
  mutate(top_genes = case_when(rownames(pbulk.markers.new) %in% rownames(top_genes) ~ "top",   
                               TRUE ~ 'ns')) 

# create a col "labels" 
pbulk.markers.new$labels <- NA
pbulk.markers.new$labels[pbulk.markers.new$top_genes != "ns"] <- pbulk.markers.new$gene[pbulk.markers.new$top_genes != "ns"]

# #----------------------------------------------------------------------------------------------
# 8. intact plot
pbulk.markers.new %>%
  ggplot(aes(x = avg_log2FC,       
             y = -log10(p_val_adj),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type,          #colour=str_wrap(gene_type,20)
             col = gene_type,
             label = labels)) + 
  geom_point(shape = 21, colour = "black") + 
  scale_color_manual(values = c("#ffad73", "#26b3ff", "#ff7387", 'grey')) +
  geom_text_repel(size = 4,  box.padding = 0.5, max.overlaps = 50) +             # box.padding = 0.5
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) +
  scale_x_continuous(breaks = c(seq(-3, 6, 2)),          # -10, 10, 2    
                     limits = c(-3, 2)) 

# 9. save in PNG/PDF
path_2 <- '/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/romans/figures/'
ggsave(paste0(path_2, "myeloid/whole/DC2/volcano_DC2_mast.pdf"),
       width = 10,    #15
       height = 6,    #10
       dpi = 300)

#--------------------------------------3. PATHWAY ANALYSIS (GSEA)----------------------------

de <- pbulk.markers 
de <- pbulk.markers %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(cell_type == 'CD36+ cDC2')
  

stats <- de$p_val_adj
names(stats) <- de$gene

load("keggSymbolHuman.rdata")       # keggSymbolHuman.rdata
fgseaResults <- fgseaMultilevel(keggSymbolHuman, stats,
                                minSize = 1, maxSize = 500,
                                nPermSimple = 10000)
# head(fgseaResults, 3)

topPathwaysUp <- fgseaResults[ES > 0, ][head(order(padj), n=5), pathway]
topPathwaysDown <- fgseaResults[ES < 0, ][head(order(pval), n=5), pathway]
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
