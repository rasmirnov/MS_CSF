setwd("~/Downloads/MD_proj/neuro_proj/data/data")
# install.packages('svglite')
library(svglite)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pbkrtest)
library(ggpubr)
library(rstatix)
library(FSA)

# 1.upload data (if duplicated - remove in scanpy)
meta <- read.csv("~/Downloads/MD_proj/neuro_proj/data/data/cd8_meta_comb.csv",
                 row.names = 1)
adata <- meta

#--------------------------------------------filtering------------------------------------------------------------

# leave only HC & MS for comparation
adata_meta <- adata %>% 
  filter(disease_group_comb == c('control', 'MS', 'ND', 'OID'))   # 'control', 'MS', 'ND', 'infectious', 'OID'

# # # # remove useless cell types 
# adata_meta <- adata_meta %>%
#   filter(!(cell_type %in% c('HLA-DR+', 'NKT-like')) | !(organ %in% c('PBMC', 'CSF')) | !(disease_group_comb %in% c('ND')))
# 
# adata_meta <- adata_meta %>%
#   filter(!(cell_type %in% c('Neutrophil')))
# #-----------------------------------------------------------------------------------------------------------

# 2.slice useless cols
adata_meta <- adata_meta %>%
  select(disease_group_comb, cell_type, organ, sample)

# 3.calculate #cells per sample
num_tot_cells <- adata_meta %>% group_by(sample) %>% 
  summarise(total_cells = n()) 

###----------------------------------------------------------------------------------------------
# 4.0
adata_meta <- adata_meta %>% 
  mutate(disease_group_comb = recode(disease_group_comb,
                                   'control' = 'HC'))

adata_meta$disease_organ <- paste0(adata_meta$disease_group_comb, '_', adata_meta$organ)

# 4.1.create list for each disease_group with unique samples: 
sampl_keep <- sapply(split(adata_meta$sample, adata_meta$disease_organ), unique)
# sampl_keep

# 4.2.count freq for all possible combination
df <- as.data.frame(table(adata_meta$disease_group_comb,
                          adata_meta$cell_type,
                          adata_meta$organ,
                          adata_meta$sample))
# df

# rename colnames
df <- rename(df, disease_group_comb = Var1, cell_type = Var2,
             organ = Var3, sample = Var4, cells_per_sample = Freq)
# df

df$disease_organ <- paste0(df$disease_group_comb, '_', df$organ)

# 4.3.exclude rows if sample was not initially presented in a disease_group
# Convert disease_group to character vector
df$disease_organ <- as.character(df$disease_organ)

# 4.4.Filter rows based on sampl_keep
cell_type_counts <- subset(df, disease_organ %in% names(sampl_keep) & 
                             unlist(Map(function(x, y) x %in% y,
                                        df$sample, sampl_keep[df$disease_organ])))

# 4.5.repeat filtering of useless cell types:
# cell_type_counts <- cell_type_counts %>%
#   filter(!(cell_type %in% c('Transitional')) | !(organ %in% c('PBMC', 'CSF')) | !(disease_group_comb %in% c('OID')))
# 
# cell_type_counts
#-----------------------------------------------------------------------------------------------

# 5.map sum column to adata_meta_per:
cell_type_counts$total_cells <- with(cell_type_counts,
                                     num_tot_cells$total_cells[match(sample, num_tot_cells$sample)])

# 6.calculate percent & rename 2 values
cell_type_counts <- cell_type_counts %>% 
  mutate(percent = cells_per_sample/total_cells * 100) %>% 
  mutate(disease_group_comb = recode(disease_group_comb,
                                     'control' = 'HC')) # 'control' = 'HC', 'infectious' = 'infec'

# 7.fix the order of categories:
cell_type_counts$disease_group_comb = factor(cell_type_counts$disease_group_comb,
                                             levels=c('HC', 'MS', 'ND', 'OID'))  # , 'HC', 'MS', 'ND','infec', 'OID'

cell_type_counts$organ = factor(cell_type_counts$organ,
                                levels=c('PBMC', 'CSF'))

cell_type_counts$cell_type = factor(cell_type_counts$cell_type,
                                    levels=c('Naïve-IFN', 'Tem GZMB+', 'Tem GZMK+',   
                                             'HLA-DR+', 'MAIT CD8+', 'NKT-like'))

# cell_type_counts <- cell_type_counts %>%
#   filter(!(cell_type %in% c('HLA-DR+', 'NKT-like')) | !(organ %in% c('PBMC', 'CSF')) | !(disease_group_comb %in% c('ND')))

# #------------------------------------------Optional---------------------------------------------------------------
# ### filter samples based on lower limit: 5%
# cell_type_counts_filtered <- cell_type_counts %>% 
#   group_by(disease_group_comb) %>%
#   filter(total_cells > quantile(total_cells, probs = 0.1))
# 
# # basic stats after filtering
# # length(unique(cell_type_counts_filtered$sample))
# # sum(cell_type_counts_filtered$cells_per_sample)
# 
# cell_type_counts <- cell_type_counts_filtered
# 
# #----------------------------------------------------------------------------
# # 7.1. create a dataframe for samples distribution
# num_tot_cells$disease_group_comb <- with(num_tot_cells,
#                                          cell_type_counts$disease_group_comb[match(sample, cell_type_counts$sample)]) #adata_meta
# num_tot_cells$organ <- with(num_tot_cells,
#                             cell_type_counts$organ[match(sample, cell_type_counts$sample)])
# # rename two groups
# num_tot_cells <- num_tot_cells %>% 
#   mutate(disease_group_comb = recode(disease_group_comb,
#                                      'control' = 'HC'))    #'control' = 'HC', 'infectious' = 'infec'
# 
# # fix the order of categories:
# num_tot_cells$disease_group_comb <-  factor(num_tot_cells$disease_group_comb,
#                                             levels=c('ND', 'OID')) # 'HC', 'MS', 'ND', 'infec'
# num_tot_cells$organ <-  factor(num_tot_cells$organ,
#                                levels=c('PBMC', 'CSF')) 
# 
# 
# # 7.2 check the distribution type:
# num_tot_cells %>% 
#   group_by(disease_group_comb, organ) %>%    # "organ" to add?
#   summarise(p_value = shapiro.test(total_cells)$p.value)
# 
# # num_tot_cells %>% 
# #   group_by(disease_group_comb) %>%    # "organ" to add?
# #   summarise(p_value = shapiro.test(total_cells)$p.value)
# 
# # 7.3 visualize it
# ggplot(num_tot_cells, aes(x=total_cells)) + 
#   geom_histogram(color="black", fill="white", alpha=0.5, bin=30) + 
#   xlab("# of cells per disease group") +
#   ylab("frequency") +
#   theme(text = element_text(size=10)) + 
#   facet_wrap(~disease_group_comb + organ, nrow=2, scales = "free")
# 
# # save
# ggsave(paste0("~/Downloads/MD_proj/neuro_proj/results/figures/histogram/.png"),
#        width = 10,
#        height = 5,
#        dpi = 300)
# 
# #--------------------------------------------------------------------------------------------------------------
# ONLY IF THE DISTRIBUTION IS HYBRIDE (NORM & NOT NORM)
# NORMAL
cell_type_counts_norm <- cell_type_counts %>%
  filter(disease_group_comb == 'ND' | disease_group_comb == 'OID')
# NOT NORMAL
cell_type_counts_unnorm <- cell_type_counts %>%
  filter(disease_group_comb == 'HC' | disease_group_comb == 'MS')

# # 7.4. tukey test (NORMAL):
tukey_stat <- cell_type_counts_norm %>%
  group_by(cell_type, disease_group_comb) %>%
  tukey_hsd(percent ~ organ) %>%
  add_xy_position(x = "cell_type") %>%
  select(disease_group_comb, cell_type, group1, group2, p.adj, p.adj.signif, y.position, groups, xmin, xmax)

# #7.6. dunn_test (NOT NORMAL):
dunn_stat <- cell_type_counts_unnorm %>%
  group_by(cell_type, disease_group_comb) %>%
  dunn_test(percent ~ organ, p.adjust.method = 'fdr') %>%
  add_xy_position(x = "cell_type") %>%            # change asterics positions here
  select(disease_group_comb, cell_type, group1, group2, p.adj, p.adj.signif, y.position, groups, xmin, xmax)

# merge two statistics together
stat <- bind_rows(dunn_stat, tukey_stat)
#---------------------------------------------------------------------------------------------------------
### REGULAR TEST
# remove organ specific populations Microglia & Neutrophil
# cell_type_counts %>%
#   filter((disease_group_comb == 'OID') & (organ == 'PBMC')) %>%
#   distinct(cell_type)
# 
# #7.6. dunn_test (NOT NORMAL):
# dunn_stat <- cell_type_counts %>%
#   group_by(cell_type, disease_group_comb) %>%
#   dunn_test(percent ~ organ, p.adjust.method = 'fdr') %>%
#   add_xy_position(x = "cell_type")          # change asterics positions here

# # 7.4. tukey test (NORMAL):
# tukey_stat <- cell_type_counts %>%
#   group_by(cell_type, disease_group_comb) %>%
#   tukey_hsd(percent ~ organ) %>%
#   add_xy_position(x = "cell_type")
# 
# tukey_stat[[11, 11]] <- 'ns'

# 8.plot 
# ggplot(cell_type_counts, aes(x = disease_group_comb, y = percent, fill = disease_group_comb))+
ggboxplot(cell_type_counts, x = 'cell_type', y = 'percent', fill = 'organ') +
  # geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),
             aes(fill = organ), colour="black", pch=21) + 
  xlab("Cell type") +
  ylab(paste0("Percentage of all CD8 T cells")) +
  stat_pvalue_manual(stat, hide.ns = T) +
  theme_grey() +
  theme(legend.position="top") + 
  guides(fill = guide_legend(nrow = 1)) +             # important to specify colour/fill
  scale_fill_discrete(name = "Organ") +
  facet_wrap(~disease_group_comb, nrow=2, scales = "free")

# 9.save
ggsave(paste0("~/Downloads/MD_proj/neuro_proj/results/figures/boxplots/cd8/PBMC_CSF_HC_MS_ND_OID.png"),
       width = 14,     # 8
       height = 6,     # 6
       dpi = 300)

# meta <- meta %>% 
#   rename(cell_type = cell_type_2, cell_type_2 = cell_type)
# 
# meta <- meta %>%
#   mutate(cell_type = replace(cell_type, cell_type == 'Naive', 
#                              'Naïve-INF'))
# write.csv(meta, 'cd4_meta.csv')      