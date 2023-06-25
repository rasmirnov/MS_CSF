setwd("~/Downloads/MD_proj/neuro_proj/data/data")
library(svglite)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pbkrtest)
library(ggpubr)
library(rstatix)
library(FSA)

# 1.upload data (if duplicated - remove in scanpy)
meta <- read.csv("~/Downloads/MD_proj/neuro_proj/data/data/atlas_meta_comb.csv",
                 row.names = 1)
adata_meta <- meta

# adata_meta <- adata_meta %>%
#   filter(cell_type %in% c('CD14 Mono', 'CD16 Mono', 'IFN CD14 Mono', 'Microglia',
#                           'cDC1', 'cDC2_2', 'pDC', 'DC-LAMP+'))
#  subset by organ
adata_meta_CSF <- adata_meta %>% 
  filter(organ == 'CSF')

# lst <- list(adata_meta_PBMC)     # adata_meta_CSF

adata <- adata_meta_CSF
organ <- adata$organ
adata_meta <- adata

# 2.slice useless cols
adata_meta <- adata_meta %>%
  select(disease_group_comb, cell_type, sample)

#-----------------------------------------------------------------------------------------------
# # count # of samples per disease group
# adata_meta %>% group_by(disease_group_comb) %>% 
#   summarise(n_samples = n_distinct(sample))
#-----------------------------------------------------------------------------------------------

# 3.calculate #cells per sample
num_tot_cells <- adata_meta %>% group_by(sample) %>% 
  summarise(total_cells = n()) 

# # 4.calculate total #cells in each disease_group_comb 
# cell_type_counts <- adata_meta %>% 
#   group_by(disease_group_comb, cell_type, sample) %>% 
#   summarise(cells_per_sample = n())    

###----------------------------------------------------------------------------------------------
# 4.1.create list for each disease_group with unique samples: 
sampl_keep <- sapply(split(adata_meta$sample, adata_meta$disease_group_comb), unique)
# sampl_keep

# 4.2.count freq for all possible combination
df <- as.data.frame(table(adata_meta$disease_group_comb,
                                        adata_meta$cell_type,
                                        adata_meta$sample))
# df

# rename colnames
df <- rename(df, disease_group_comb = Var1, cell_type = Var2,
             sample = Var3, cells_per_sample = Freq)
# df

# 4.3.exclude rows if sample was not initially presented in a disease_group
# Convert disease_group to character vector
df$disease_group_comb <- as.character(df$disease_group_comb)

# 4.4.Filter rows based on sampl_keep
cell_type_counts <- subset(df, disease_group_comb %in% names(sampl_keep) & 
                        unlist(Map(function(x, y) x %in% y,
                                   df$sample, sampl_keep[df$disease_group_comb])))
# cell_type_counts
#-----------------------------------------------------------------------------------------------

# 5.map sum column to adata_meta_per:
cell_type_counts$total_cells <- with(cell_type_counts,
                                     num_tot_cells$total_cells[match(sample, num_tot_cells$sample)])

# 6.calculate percent & rename 2 values
cell_type_counts <- cell_type_counts %>% 
  mutate(percent = cells_per_sample/total_cells * 100) %>% 
  mutate(disease_group_comb = recode(disease_group_comb, 'control' = 'HC', 'infectious' = 'infec'))

# 7.fix the order of categories:
cell_type_counts$disease_group_comb = factor(cell_type_counts$disease_group_comb,
                                             levels=c('HC', 'MS', 'ND', 'infec', 'OID'))  # 'HC', 'MS', 'ND', 'infec', 'OID'

cell_type_counts$cell_type = factor(cell_type_counts$cell_type,
                                    levels=c('Myeloid cells', 'CD4 T cells', 'CD8 T cells',  
                                             'B cells', 'NK cells'))

# #------------------------------------------Optional---------------------------------------------------------------
# ### filter samples based on lower limit: 5%
# cell_type_counts_filtered <- cell_type_counts %>% 
#   group_by(disease_group_comb) %>%
#   filter(total_cells > quantile(total_cells, probs = 0.1))
# 
# cell_type_counts <- cell_type_counts_filtered

#----------------------------------------------------------------------------
# # 7.1. create a dataframe for samples distribution
# num_tot_cells$disease_group_comb <- with(num_tot_cells,
#                                          cell_type_counts$disease_group_comb[match(sample, cell_type_counts$sample)]) #adata_meta
# # rename two groups
# num_tot_cells <- num_tot_cells %>%
#   mutate(disease_group_comb = recode(disease_group_comb, 'control' = 'HC', 'infectious' = 'infec'))
# 
# # fix the order of categories:
# num_tot_cells$disease_group_comb <-  factor(num_tot_cells$disease_group_comb,
#                                             levels=c('HC', 'MS', 'ND', 'infec', 'OID'))
# 
# # # filter out the group, if >2 sample
# # num_tot_cells <- num_tot_cells %>%
# #   filter(disease_group_comb != c('OID'))          # c('ND', 'infec')
# # num_tot_cells <- num_tot_cells %>%
# #   filter(disease_group_comb %in% c('HC', 'CIS', 'MS', 'OID'))
# 
# # 7.2 check the distribution type:
# num_tot_cells %>%
#   group_by(disease_group_comb) %>%
#   summarise(p_value = shapiro.test(total_cells)$p.value)
# # shapiro.test(cell_type_counts$total_cells)$p.value
# 
# # 7.3 visualize it
# ggplot(num_tot_cells, aes(x=total_cells)) +
#   geom_histogram(color="black", fill="white", alpha=0.5, bin=30) +
#   xlab("# of cells per disease group") +
#   ylab("frequency") +
#   theme(text = element_text(size=10)) +
#   facet_wrap(~disease_group_comb, nrow=2, scales = "free")
# 
# # save
# ggsave(paste0("~/Downloads/MD_proj/neuro_proj/results/figures/histogram/atlas", organ, ".png"),
#        width = 10,
#        height = 5,
#        dpi = 300)
# #--------------------------------------------------------------------------------------------------------------
# # ONLY IF THE DISTRIBUTION IS HYBRIDE (NORM & NOT NORM)
# # NORMAL
# cell_type_counts_norm <- cell_type_counts %>%
#   filter(disease_group_comb == 'HC' | disease_group_comb == 'infec' | disease_group_comb == 'OID')
# # NOT NORMAL
# cell_type_counts_unnorm <- cell_type_counts %>%
#   filter(disease_group_comb == 'HC' | disease_group_comb == 'MS' | disease_group_comb == 'ND')
# 
# # remove neutrophils
# # 7.4. tukey test (NORMAL):
# tukey_stat <- cell_type_counts_norm %>%
#   group_by(cell_type) %>%
#   tukey_hsd(percent ~ disease_group_comb) %>%
#   filter(group1 == 'HC') %>%
#   add_xy_position(x = "disease_group_comb") %>%
#   select(cell_type, group1, group2, p.adj, p.adj.signif, y.position, groups, xmin, xmax)                        # leave only specific colms
# 
# # 7.5. dunn_test (NOT NORMAL):
# dunn_stat <- cell_type_counts_unnorm %>%
#   group_by(cell_type) %>%
#   dunn_test(percent ~ disease_group_comb, p.adjust.method = 'fdr') %>%  # give an error if >2 groups to compare
#   filter(group1 == c('HC')) %>%
#   add_xy_position(x = "disease_group_comb") %>%
#   select(cell_type, group1, group2, p.adj, p.adj.signif, y.position, groups, xmin, xmax)
# 
# # optional: change the position of sign.усиков
# dunn_stat[[15, 6]] <- 20
# 
# # merge two statistics together
# stat <- bind_rows(dunn_stat, tukey_stat)
#---------------------------------------------------------------------------------------------------------
### REGULAR TEST
# # 7.6. dunn_test (NOT NORMAL):
dunn_stat <- cell_type_counts %>%
  group_by(cell_type) %>%
  dunn_test(percent ~ disease_group_comb, p.adjust.method = 'fdr') %>%
  filter(group1 == c('HC')) %>%
  add_xy_position(x = "disease_group_comb")

# # 7.4. tukey test (NORMAL):
# tukey_stat <- cell_type_counts %>%
#   group_by(cell_type) %>%
#   tukey_hsd(percent ~ disease_group_comb) %>%
#   filter(group1 == 'HC') %>%
#   add_xy_position(x = "disease_group_comb") %>%
#   select(cell_type, group1, group2, p.adj, p.adj.signif, y.position, groups, xmin, xmax)

# 8.plot 
# ggplot(cell_type_counts, aes(x = disease_group_comb, y = percent, fill = disease_group_comb))+
ggboxplot(cell_type_counts, x = 'disease_group_comb', y = 'percent', fill = 'disease_group_comb') +
  # geom_boxplot() +
  geom_point(aes(fill = disease_group_comb), colour="black", pch=21) + 
  xlab("Disease group") +
  ylab(paste0("Percentage of all ", organ,  " cells")) +
  # stat_compare_means(method = 'kruskal.test', label.x = 2.3, size = 2.8) + 
  stat_pvalue_manual(dunn_stat, hide.ns = T) +
  theme_grey() +
  theme(legend.position="top") + 
  guides(fill = guide_legend(nrow = 1)) +             # important to specify colour/fill
  scale_fill_discrete(name = "Disease group") +
  facet_wrap(~cell_type, nrow=2, scales = "free")     # nrow=2

# 9.save
ggsave(paste0("~/Downloads/MD_proj/neuro_proj/results/figures/boxplots/atlas/", organ, ".png"),
       width = 9,    #9 or 12
       height = 5.5,    #6
       dpi = 300)


# to check color names 
# ggplot_build(p)$data
# HC: #F8766D
# MS: #A3A500
# ND: 00B0F6
# infec: 00BF7D
# OID: E76BF3


# meta <- meta %>% 
#   rename(cell_type = cell_type_2, cell_type_2 = cell_type)
# 
# meta <- meta %>%
#   mutate(disease_group_comb = replace(disease_group, disease_group == 'CIS',
#                              'MS'))
# write.csv(meta, 'myeloid_meta_comb.csv')
