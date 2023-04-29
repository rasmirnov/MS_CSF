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
adata <- meta

# 2.slice useless cols
adata_meta <- adata_meta %>%
  select(cell_type, organ, sample)

# 3.calculate #cells per sample
num_tot_cells <- adata_meta %>% group_by(sample) %>% 
  summarise(total_cells = n()) 

#----------------------------------------------------------------------------
# 4.calculate total #cells in each disease_group_comb 

# 4.1.create list for each disease_group with unique samples: 
sampl_keep <- sapply(split(adata_meta$sample, adata_meta$organ), unique)
# sampl_keep

# 4.2.count freq for all possible combination
df <- as.data.frame(table(adata_meta$organ,
                          adata_meta$cell_type,
                          adata_meta$sample))
# df

# rename colnames
df <- rename(df, organ = Var1, cell_type = Var2,
             sample = Var3, cells_per_sample = Freq)
# df

# 4.3.exclude rows if sample was not initially presented in a disease_group
# Convert disease_group to character vector
df$organ <- as.character(df$organ)

# Filter rows based on sampl_keep
cell_type_counts <- subset(df, organ %in% names(sampl_keep) & 
                             unlist(Map(function(x, y) x %in% y,
                                        df$sample, sampl_keep[df$organ])))
# cell_type_counts

#----------------------------------------------------------------------------

# 5.map sum column to adata_meta_per:
cell_type_counts$total_cells <- with(cell_type_counts,
                                     num_tot_cells$total_cells[match(sample, num_tot_cells$sample)])

# 6.calculate percent & rename 2 values
cell_type_counts <- cell_type_counts %>% 
  mutate(percent = cells_per_sample/total_cells * 100)

# 7.fix the order of categories:

cell_type_counts$organ = factor(cell_type_counts$organ,
                                levels=c('PBMC', 'CSF'))

cell_type_counts$cell_type = factor(cell_type_counts$cell_type,
                                    levels=c('Myeloid cells', 'CD4 T cells', 'CD8 T cells',  
                                             'B cells', 'NK cells'))

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
# num_tot_cells$organ <- with(num_tot_cells,
#                                          cell_type_counts$organ[match(sample, cell_type_counts$sample)]) #adata_meta
# 
# # fix the order of categories:
# num_tot_cells$organ <-  factor(num_tot_cells$organ,
#                                levels=c('PBMC', 'CSF'))
# 
# # 7.2 check the distribution type:
# num_tot_cells %>%
#   group_by(organ) %>%    # "organ" to add?
#   summarise(p_value = shapiro.test(total_cells)$p.value)
# 
# # num_tot_cells %>%
# #   group_by(disease_group_comb) %>%    # "organ" to add?
# #   summarise(p_value = shapiro.test(total_cells)$p.value)
# 
# # 7.3 visualize it
# ggplot(num_tot_cells, aes(x=total_cells)) +
#   geom_histogram(color="black", fill="white", alpha=0.5, bin=30) +
#   xlab("# of cells") +
#   ylab("frequency") +
#   theme(text = element_text(size=10)) +
#   facet_wrap(~organ, nrow=2, scales = "free")
# 
# # save
# ggsave(paste0("~/Downloads/MD_proj/neuro_proj/results/figures/histogram/.png"),
#        width = 10,
#        height = 5,
#        dpi = 300)

# #--------------------------------------------------------------------------------------------------------------
# # ONLY IF THE DISTRIBUTION IS HYBRIDE (NORM & NOT NORM)
# # NORMAL
# cell_type_counts_norm <- cell_type_counts %>%
#   filter(disease_group_comb == 'MS' | disease_group_comb == 'OID')
# # NOT NORMAL
# cell_type_counts_unnorm <- cell_type_counts %>%
#   filter(disease_group_comb == 'HC' | disease_group_comb == 'ND')
# 
# # # 7.4. tukey test (NORMAL):
# tukey_stat <- cell_type_counts_norm %>%
#   group_by(cell_type, disease_group_comb) %>%
#   tukey_hsd(percent ~ organ) %>%
#   add_xy_position(x = "cell_type") %>%
#   select(disease_group_comb, cell_type, group1, group2, p.adj, p.adj.signif, y.position, groups, xmin, xmax)
# 
# # #7.6. dunn_test (NOT NORMAL):
# dunn_stat <- cell_type_counts_unnorm %>%
#   group_by(cell_type, disease_group_comb) %>%
#   dunn_test(percent ~ organ, p.adjust.method = 'fdr') %>%
#   add_xy_position(x = "cell_type") %>%            # change asterics positions here
#   select(disease_group_comb, cell_type, group1, group2, p.adj, p.adj.signif, y.position, groups, xmin, xmax)
# 
# # merge two statistics together
# stat <- bind_rows(dunn_stat, tukey_stat)
# # #---------------------------------------------------------------------------------------------------------
# ### REGULAR TEST
# # remove organ specific populations Microglia & Neutrophil
# cell_type_counts %>%
#   filter((disease_group_comb == 'ND') & (organ == 'CSF')) %>%
#   distinct(cell_type)
# 
#7.6. dunn_test (NOT NORMAL):
dunn_stat <- cell_type_counts %>%
  group_by(cell_type) %>%
  dunn_test(percent ~ organ, p.adjust.method = 'fdr') %>%
  add_xy_position(x = "cell_type")          # change asterics positions here

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
  ylab(paste0("Percentage of all cells")) +
  stat_pvalue_manual(dunn_stat, hide.ns = T) +
  theme_grey() +
  theme(legend.position="top") + 
  guides(fill = guide_legend(nrow = 1)) +             # important to specify colour/fill
  scale_fill_discrete(name = "Organ") 

# 9.save
ggsave(paste0("~/Downloads/MD_proj/neuro_proj/results/figures/boxplots/atlas/PBMC_CSF_all.svg"),
       width = 7,     # 8
       height = 4,     # 6
       dpi = 300)
