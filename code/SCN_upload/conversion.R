# in plasma_cells folder now
install.packages("stringi", dependencies = TRUE, repos = "http://cran.us.r-project.org")
library(stringi)
library(SCNPrep) 
library(Seurat)
library(data.table)
library(Matrix)
library(dplyr)
library(readr)

clusters <- c('cd4_t_cells')

for (cluster in clusters) {
  
  setwd(paste0('/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/romans/csf_atlas_new/',
               cluster, '/no_duplicates_new'))
  object <- readRDS(paste0(cluster, '.rds'))
  sample_id <- cluster
  tokenn <- paste(c(sample_id, stri_rand_strings(1, 8)), collapse="_")
  # dir.create('plasma_cells')
  
  markers_files <- list.files(pattern = 'markers.tsv', recursive = T, full.names = T)
  markers <- lapply(markers_files, fread)
  names(markers) <- gsub('/|markers', '', stringr::str_extract(markers_files, 'markers/*.*/'))
  # change colnames: avg_log2FC to avg_logFC: list with new names
  markers_upd <- lapply(markers, function(i) {gsub('avg_log2FC', 'avg_logFC', colnames(i))})
  # replace old_names on new_names:
  for (i in seq_along(markers)){
    colnames(markers[[i]]) <- markers_upd[[i]]
  }
  
  object@meta.data <- object@meta.data %>%
    mutate_if(is.character, as.factor)  
  
  migrateSeuratObject(object,
                      assay="RNA",
                      species='hs',
                      outdir = tokenn,
                      public = F,
                      curated = F,
                      markers = markers,
                      generateMarkers = F,
                      name= sample_id,  
                      token= tokenn)
}

