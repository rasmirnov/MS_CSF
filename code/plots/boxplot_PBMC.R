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

#--------------------------------------------filtering------------------------------------------------------------
#  subset by organ
adata_meta_PBMC <- adata_meta %>% 
  filter(organ %in% 'PBMC')

# leave only HC & MS for comparation
adata_meta_PBMC <- adata_meta_PBMC %>% 
  filter(disease_group_comb %in% c('control', 'MS'))

# # remove Neutrophil (HC) & DC-LAMP+ (MS) because they're present in only one group
# adata_meta_PBMC <- adata_meta_PBMC %>%
#   filter(cell_type %in% c('CD14 Mono', 'IFN CD14 Mono', 'CD16 Mono',
#                           'cDC1', 'cDC2_2', 'pDC', 'DC-LAMP+'))
# #-----------------------------------------------------------------------------------------------------------

adata_meta <- adata_meta_PBMC
organ <- adata_meta$organ

# 2.slice useless cols
adata_meta <- adata_meta %>%
  select(disease_group_comb, cell_type, sample)

# 3.calculate #cells per sample
num_tot_cells <- adata_meta %>% group_by(sample) %>% 
  summarise(total_cells = n()) 

#----------------------------------------------------------------------------
# 4.calculate total #cells in each disease_group_comb 

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

# Filter rows based on sampl_keep
cell_type_counts <- subset(df, disease_group_comb %in% names(sampl_keep) & 
                             unlist(Map(function(x, y) x %in% y,
                                        df$sample, sampl_keep[df$disease_group_comb])))
# cell_type_counts

#----------------------------------------------------------------------------

# 5.map sum column to adata_meta_per:
cell_type_counts$total_cells <- with(cell_type_counts,
                                     num_tot_cells$total_cells[match(sample, num_tot_cells$sample)])

# 6.calculate percent & rename 2 values
cell_type_counts <- cell_type_counts %>%
  mutate(percent = cells_per_sample/total_cells * 100) %>%
  mutate(disease_group_comb = recode(disease_group_comb, 'control' = 'HC'))

# 7.fix the order of categories:
cell_type_counts$disease_group_comb = factor(cell_type_counts$disease_group_comb,
                                             levels=c('HC', 'MS'))  # , 'CIS', 'infec', 'ND', 'OID'

cell_type_counts$cell_type = factor(cell_type_counts$cell_type,
                                    levels=c('Myeloid cells', 'CD4 T cells', 'CD8 T cells',  
                                             'B cells', 'NK cells'))

# #------------------------------------------Optional---------------------------------------------------------------
# ### filter samples based on lower limit: 5%
# cell_type_counts_filtered <- cell_type_counts %>% 
#   group_by(disease_group_comb) %>%
#   filter(total_cells > quantile(total_cells, probs = 0.2))
# 
# # basic stats after filtering
# # length(unique(cell_type_counts_filtered$sample))
# # sum(cell_type_counts_filtered$cells_per_sample)
# 
# cell_type_counts <- cell_type_counts_filtered

#----------------------------------------------------------------------------
# # 7.1. create a dataframe for samples distribution
# num_tot_cells$disease_group_comb <- with(num_tot_cells,
#                                          cell_type_counts$disease_group_comb[match(sample, cell_type_counts$sample)]) #adata_meta
# # rename two groups
# num_tot_cells <- num_tot_cells %>%
#   mutate(disease_group_comb = recode(disease_group_comb, 'control' = 'HC'))
# 
# # fix the order of categories:
# num_tot_cells$disease_group_comb <-  factor(num_tot_cells$disease_group_comb,
#                                             levels=c('HC', 'MS'))
# 
# 
# # 7.2 check the distribution type:
# num_tot_cells %>%
#   group_by(disease_group_comb) %>%
#   summarise(p_value = shapiro.test(total_cells)$p.value)
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
# ggsave(paste0("~/Downloads/MD_proj/neuro_proj/results/figures/histogram/atlas_", organ, ".png"),
#        width = 10,
#        height = 5,
#        dpi = 300)
# # #---------------------------------------------------------------------------------------------------------
# ### REGULAR TEST
# 
# cell_type_counts %>%
#   filter(disease_group_comb == 'HC') %>%
#   distinct(cell_type)

# 7.6. dunn_test (NOT NORMAL):
dunn_stat <- cell_type_counts %>%
  group_by(cell_type) %>%
  dunn_test(percent ~ disease_group_comb, p.adjust.method = 'fdr') %>%
  filter(group1 == c('HC')) %>%
  add_xy_position(x = "cell_type")

# # # 7.4. tukey test (NORMAL):
# tukey_stat <- cell_type_counts %>%
#   group_by(cell_type) %>%
#   tukey_hsd(percent ~ disease_group_comb) %>%
#   filter(group1 == 'HC') %>%
#   add_xy_position(x = "cell_type")

# 8.plot (change the color for MS on #A3A500)
# ggplot(cell_type_counts, aes(x = disease_group_comb, y = percent, fill = disease_group_comb))+
ggboxplot(cell_type_counts, x = 'cell_type', y = 'percent', fill = 'disease_group_comb') +
  # geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),
             aes(fill = disease_group_comb), colour="black", pch=21) + 
  xlab("Cell type") +
  ylab(paste0("Percentage of all ", organ,  " cells")) +
  stat_pvalue_manual(dunn_stat, hide.ns = T) +
  theme_grey() +
  theme(legend.position="top") + 
  guides(fill = guide_legend(nrow = 1)) +             # important to specify colour/fill
  scale_fill_manual(name = "Disease group",
                    values = c("HC" = "#F8766D",
                               "MS" = "#A3A500"))

# 9.save
ggsave(paste0("~/Downloads/MD_proj/neuro_proj/results/figures/boxplots/atlas/", organ, ".png"),
       width = 8,     # 10
       height = 4,     # 5
       dpi = 300)
