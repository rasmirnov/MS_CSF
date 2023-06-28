library(Seurat)
library(dplyr)
library(ggplot2)

# upload data (if duplicated - remove in scanpy)
atlas_meta <- read.csv("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/romans/csf_atlas_new/myeloid_cells/whole/myeloid_meta_comb.csv",
                       row.names = 1)
adata_meta <- atlas_meta

# slice useless cols
adata_meta <- adata_meta %>%
  select(cell_type, disease_group, organ) %>% 
  filter(cell_type %in% c('CD36+ cDC2', 'CD32B+ cDC2 like'))   ### optional

# split by organ
adata_meta_CD36 <- adata_meta[adata_meta$cell_type == 'CD36+ cDC2', ]
adata_meta_CD32 <- adata_meta[adata_meta$cell_type == 'CD32B+ cDC2 like', ]

# calculate percentages: PBMC
adata_meta_CD36 <- adata_meta_CD36 %>% group_by(organ) %>%  #disease_group,
  summarise(per_cluster = n()) %>%
  mutate(sum_cells = sum(per_cluster)) %>%
  mutate(percent = per_cluster/sum_cells*100) %>% 
  mutate(cell_type = 'CD36+ cDC2') %>% 
  mutate(across(where(is.numeric), round, 1))
  # mutate(organ = 'PBMC')

# calculate percentages: CSF
adata_meta_CD32 <- adata_meta_CD32 %>% group_by(organ, cell_type) %>%    #disease_group, 
  summarise(per_cluster = n()) %>%
  mutate(sum_cells = sum(per_cluster)) %>%
  mutate(percent = per_cluster/sum_cells*100)  
  # mutate(organ = 'CSF')

# merge them
adata_meta_per <- rbind(adata_meta_CD36, adata_meta_CD32)

# fix the order of categories:
# adata_meta_per$disease_group = factor(adata_meta_per$disease_group, levels=c('control', 'CIS', 'MS', 'ND', 'infectious', 'OID'))
adata_meta_per$cell_type = factor(adata_meta_per$cell_type, levels=c('CD36+ cDC2', 'CD32B+ cDC2 like'))
adata_meta_per$organ = factor(adata_meta_per$organ, levels=c('PBMC', 'CSF'))

# plot
ggplot(adata_meta_per, aes(x = cell_type, y = percent, fill = as.factor(organ), label = percent))+     # disease_group
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "Organ") +
  geom_text(size = 4.5, position = position_stack(vjust = 0.5)) +
  # xlab("Cell type") +
  ylab("Percentage of cells") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        axis.title.x = element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=13))
  # facet_grid(rows = vars(organ))

# save
ggsave("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/romans/figures/myeloid/whole/DC2_stacked_barplot.pdf",
       width = 5,
       height = 5,
       dpi = 300)

