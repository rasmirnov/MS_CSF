suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(magrittr))
suppressMessages(library(sctransform))
suppressMessages(library(Seurat))
suppressMessages(library(argparse))

set.seed(1)


parser <-
  ArgumentParser(description = 'Get scRNA-seq related figures from the paper')
parser$add_argument('--meta',
                    type = "character",
                    help = 'Path to SRA metadata anntotation')
parser$add_argument('--out_dir',
                    type = "character",
                    help = 'Path to output directory')

## SET VARIABLES

args <- parser$parse_args()

print(args)

## FUNCTIONS

add_metadata <- function(data) {
  mito.genes <-
    grep(pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-",
         x = rownames(x = GetAssayData(object = data)),
         value = TRUE)
  percent.mito <-
    Matrix::colSums(GetAssayData(object = data, slot = "counts")[mito.genes, ]) /
    Matrix::colSums(GetAssayData(object = data, slot = "counts"))
  data[['percent.mito']] <- percent.mito
  data[['percent.mito_log10']] <- log10(data[['percent.mito']] + 1)
  data[['nCount_RNA_log10']] <- log10(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log10']] <- log10(data[['nFeature_RNA']] + 1)
  data[['nCount_RNA_log2']] <- log2(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log2']] <- log2(data[['nFeature_RNA']] + 1)
  data[['scaled_mito']] <- scale(percent.mito)
  data[['scaled_nCount_RNA']] <- scale(data[['nCount_RNA_log10']])
  attr(data$scaled_nCount_RNA, "scaled:center") <- NULL
  attr(data$scaled_nCount_RNA, "scaled:scale") <- NULL
  attr(data$scaled_mito, "scaled:center") <- NULL
  attr(data$scaled_mito, "scaled:scale") <- NULL
  data
}

get_conf_interval <- function(dataset, parameter) {
  left <- mean(dataset[[parameter]][[1]]) - qnorm(0.975)
  right <- mean(dataset[[parameter]][[1]]) + qnorm(0.975)
  return(c(left, right))
}

filter_mito <- function(dataset, path){
  mt_dist <- as.data.frame(dataset[['scaled_mito']][[1]])
  mt_pers <- as.data.frame(dataset[['percent.mito']][[1]])
  scaled_mito_percentage <- scale(dataset[['percent.mito']][[1]])
  colnames(mt_dist) <- 'scaled_mito'
  colnames(mt_pers) <- 'percent.mito'
  filtration_coord <- get_conf_interval(dataset, 'scaled_mito')[2] * 
    attr(scaled_mito_percentage, 'scaled:scale') + 
    attr(scaled_mito_percentage, 'scaled:center')
  expr <- FetchData(object = dataset, vars = 'scaled_mito')
  dataset <- dataset[, which(x = expr < get_conf_interval(dataset, 'scaled_mito')[2])]
  dataset
}

get_df <- function(path) {
  data <- Read10X(path)
  data <- CreateSeuratObject(data, min.cells = 3, min.features = 200)
  data <- add_metadata(data)
  folds <- sapply(meta$folder, function(x) stringr::str_extract(path, x))
  folds <- folds[!is.na(folds)]
  data$sample <- meta$sample[match(folds, meta$folder)]
  data$disease <- meta$disease[match(folds, meta$folder)]
  data$tissue <- meta$tissue[match(folds, meta$folder)]
  data$sex <- meta$sex[match(folds, meta$folder)]
  data$age <- meta$age[match(folds, meta$folder)]
  data
}

get_data <- function(pathes) {
  objects <- lapply(pathes, get_df)
  names(objects) <- sapply(objects, function(x) unique(x$sample))
  objects
}


## GATHERING DATA TOGETHER


options(future.globals.maxSize = 15000 * 1024^2)

meta <- fread(args$meta)

pathes <- gsub('barcodes.tsv.gz', '', list.files('./', pattern = 'barcodes.tsv.gz', full.names = T, recursive = T))
print(pathes)

whole <- get_data(pathes)
setwd(args$out_dir)

## Number of cells before

cells.before <- sapply(whole, function(x) dim(GetAssayData(object = x, slot = "counts"))[2])

## NORMALIZATION

whole <- sapply(whole, function(x) SCTransform(
 x,
 ncells=min(100000, ncol(x)),
 vars.to.regress = c("percent.mito"),
 verbose = T,
 conserve.memory = T
))

## CELL CYCLE SCORING

s.genes <- stringr::str_to_title(cc.genes$s.genes)
g2m.genes <- stringr::str_to_title(cc.genes$g2m.genes)


whole <- sapply(whole, function(x) CellCycleScoring(x, s.features = s.genes,
                                                   g2m.features = g2m.genes, set.ident = F,
                                                   assay = 'SCT'))
## INTEGRATION

whole.features <- SelectIntegrationFeatures(object.list = whole, nfeatures = 2000)

whole <- PrepSCTIntegration(object.list = whole, anchor.features = whole.features,
                            verbose = FALSE)
whole.anchors <- FindIntegrationAnchors(object.list = whole, normalization.method = "SCT",
                                        anchor.features = whole.features, verbose = FALSE)
whole.integrated <- IntegrateData(anchorset = whole.anchors, normalization.method = "SCT",
                                  verbose = FALSE)

## PCA
gc()

whole.integrated <- RunPCA(whole.integrated, verbose = FALSE)

## UMAP

whole.integrated <- RunUMAP(whole.integrated, dims = 1:30)


## CLUSTERING

whole.integrated <- FindNeighbors(object = whole.integrated, dims = 1:30)
whole.integrated <- FindClusters(object = whole.integrated, resolution = )     # seq(0.2, 1, 0.2)

## SAVING: DATASET

save(list = c('whole.integrated', 'whole.features', 'whole.anchors'), file = "object.RData")

## FINDING ANS SAVING MARKERS

analyze_object <- function(object, ident) {
  Idents(object) <- object[[ident]]
  if (length(levels(object)) == 1) {
    return(message(sprintf('%s: since only one cluster was identified, markers can not be found', ident)))
  }
  out_dir <- paste0('markers/', ident)
  dir.create(out_dir, recursive = T)
  whole.markers <- FindAllMarkers(object = object,
                                  assay='SCT',
                                  only.pos = TRUE,
                                  min.pct = 0.10,
                                  test.use = 'wilcox',
                                  max.cells.per.ident = 3e3,
                                  random.seed = 42)
  write.table(whole.markers, paste(out_dir, "markers.tsv", sep = '/'), sep="\t", quote=F, row.names=F)
}

idents <- grep('integrated', colnames(whole.integrated@meta.data), value = T)
sapply(idents, function(ident) analyze_object(object = whole.integrated, ident = ident))


## Number of cells after

cells.after <- sapply(whole, function(x) length(colnames(x = x)))
cells.diff <- cells.before-cells.after
rbind(cells.before, cells.after, cells.diff)
