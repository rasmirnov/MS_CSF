setwd("~/Downloads/MD_proj/neuro_proj/data/data")
library(dplyr)
library(ggplot2)
library(tidyr)
library(pbkrtest)
library(ggpubr)
library(rstatix)
library(FSA)
library(gridExtra)

# 1.upload data (if duplicated - remove in scanpy)
meta <- read.csv("~/Downloads/MD_proj/neuro_proj/data/data/atlas_meta_comb.csv",
                 row.names = 1)
adata_meta <- meta

# replace some 'myeloid cells' on 'Microglia' based on "cell_type_2" colmn
adata_meta <- adata_meta %>% 
  mutate(cell_type = ifelse(cell_type == 'Myeloid cells' & 
                              cell_type_2 == 'Microglia', 'Microglia', cell_type))

adata_meta <- adata_meta %>% 
  filter(cell_type  %in% c('CD4 T cells', 'Microglia'))
#-----------------------------------------------------------------------------------------------
# # count # of samples per disease group
# cell_type_counts %>% group_by(disease_group_comb, cell_type) %>%
#   summarise(n_samples = n_distinct(sample))
#-----------------------------------------------------------------------------------------------

#  subset by organ
adata_meta_CSF <- adata_meta %>% 
  filter(organ %in% 'CSF')

adata <- adata_meta_CSF
organ <- adata$organ
adata_meta <- adata

# 2.slice useless cols
adata_meta <- adata_meta %>%
  select(disease_group_comb, cell_type, sample)

# 3.calculate #cells per sample
num_tot_cells <- adata_meta %>% group_by(sample) %>% 
  summarise(total_cells = n()) 

###----------------------------------------------------------------------------------------------
# 4.0
adata_meta <- adata_meta %>% 
  mutate(disease_group_comb = recode(disease_group_comb,
                                     'control' = 'HC', 'infectious' = 'infec'))

adata_meta$disease_cell_type <- paste0(adata_meta$disease_group_comb, '_', adata_meta$cell_type)

# 4.1.create list for each disease_group with unique samples: 
sampl_keep <- sapply(split(adata_meta$sample, adata_meta$disease_cell_type), unique)

# 4.2.count freq for all possible combination
df <- as.data.frame(table(adata_meta$disease_group_comb,
                          adata_meta$cell_type,
                          adata_meta$sample))

# rename colnames
df <- rename(df, disease_group_comb = Var1, cell_type = Var2,
             sample = Var3, cells_per_sample = Freq)

df$disease_cell_type <- paste0(df$disease_group_comb, '_', df$cell_type)

# 4.3.exclude rows if sample was not initially presented in a disease_group
df$disease_cell_type <- as.character(df$disease_cell_type)

# 4.4.Filter rows based on sampl_keep
cell_type_counts <- subset(df, disease_cell_type %in% names(sampl_keep) & 
                             unlist(Map(function(x, y) x %in% y,
                                        df$sample, sampl_keep[df$disease_cell_type])))

# 4.5.leave only duplicated samples (those which are present in both cell types):
cell_type_counts <- cell_type_counts %>% 
  group_by(sample, disease_group_comb) %>% filter(n() > 1)
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

# Bc/MG ratio of medians of a % in CSF:
cell_type_counts_flt <- cell_type_counts %>% 
  filter(cell_type == c('CD4 T cells', 'Microglia')) %>%
  group_by(disease_group_comb, sample) %>% 
  summarise(ratio_of_percent = percent[cell_type == "CD4 T cells"] / percent[cell_type == "Microglia"]) 

# plot
ggplot(cell_type_counts_flt, aes(disease_group_comb, ratio_of_percent, fill = disease_group_comb)) +
  geom_bar(position = 'dodge', stat = 'summary', fun.y = 'mean') +
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0.1) +
  geom_point(aes(x = disease_group_comb), shape = 21, position = position_dodge(width = 1)) + 
  xlab("Disease group") +
  ylab("CD4/Micro ratio of percentages in CSF") + 
  scale_fill_discrete(name = "Disease group") +
  theme(legend.title=element_text(size=9),
        axis.title=element_text(size=9))


# save
ggsave(paste0("~/Downloads/MD_proj/neuro_proj/results/figures/boxplots/atlas/cd4_MG_ratio.png"),
       width = 5,    #9 or 12
       height = 4,    #6
       dpi = 300)






