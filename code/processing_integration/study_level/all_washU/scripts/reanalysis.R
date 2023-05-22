suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(magrittr))
suppressMessages(library(sctransform))
suppressMessages(library(Seurat))
suppressMessages(library(argparse))
suppressMessages(library(glmGamPoi))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(miQC))
suppressMessages(library(flexmix))
suppressMessages(library(SCNPrep))
suppressMessages(library(DoubletFinder))
suppressMessages(library(jsonlite))


set.seed(1)


parser <-
  ArgumentParser(description = 'Get scRNA-seq related figures from the paper')
parser$add_argument('--out_dir',
                    type = "character",
                    help = 'Path to output directory')
parser$add_argument('--js_f',
                    type = "character",
                    help = 'path to json with given filters for data')
parser$add_argument('--sample_id',
                    type = "character",
                    help = 'Sample ID')
## SET VARIABLES

args <- parser$parse_args()

print(args)

## GATHERING DATA TOGETHER

options(future.globals.maxSize = 20000 * 1024^2)

add_metadata <- function(data) {
  mito.genes <-
    grep(pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-",
         x = rownames(x = GetAssayData(object = data)),
         value = TRUE)
  percent.mito <-
    Matrix::colSums(GetAssayData(object = data, slot = "counts")[mito.genes, ]) /
    Matrix::colSums(GetAssayData(object = data, slot = "counts"))
  data[['percent.mito']] <- percent.mito
  data[['percent.mt']] <- percent.mito
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

get_data <- function(data, js_f, out_dir) {
  load(data)
  meta <- whole.integrated@meta.data
  if ("target_clusters" %in% names(js_f)) {
    meta <- meta %>% filter(!!as.symbol(names(js_f[['target_clusters']])) %in%
                              strsplit(js_f[['target_clusters']][1][[1]], " ")[[1]])
  }
  whole.integrated <- subset(whole.integrated, cells = rownames(meta))
  DefaultAssay(whole.integrated) <- 'RNA'
  whole.integrated[['integrated']] <- NULL
  whole.integrated[['SCT']] <- NULL
  SplitObject(whole.integrated, split.by = 'sample')
}


get_df <- function(path) {
  data <- Read10X(path)
  df <- annot %>% filter(sample %in% stringr::str_split(path, '/')[[1]])
  obj <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
  obj <- add_metadata(obj)
  obj$study <- df$study
  obj$Disease <- df$Disease
  obj$Tissue <- df$Tissue
  obj$Sex <- df$Sex
  obj$Age <- df$Age
  obj$protocol <- df$protocol
  obj$sample <- df$sample
  obj
}

get_whole_obj <- function(pathes) {
  objects <- lapply(pathes, get_df)
  names(objects) <- sapply(objects, function(x) unique(x$sample))
  objects
}

get_seurat <- function(js_f) {
  if ("data" %in% names(js_f)) {
    pathes <- gsub('barcodes.*', '',
                   list.files(js_f$data, pattern = 'barcodes.tsv.gz',
                              full.names = T, recursive = T))
    unprocessed_smpls <- get_whole_obj(pathes)
  }
  if ("object" %in% names(js_f)) {
    processed_smpls <- get_data(data = js_f$object, js_f = js_f,
                                out_dir = args$out_dir)
  }
  if ("data" %in% names(js_f) & "object" %in% names(js_f)) {
    return(c(processed_smpls, unprocessed_smpls))
  }
  if ("data" %in% names(js_f)) {
    return(unprocessed_smpls)
  }
  if ("object" %in% names(js_f)) {
    return(processed_smpls)
  }
}

js_f <- jsonlite::read_json(args$js_f, simplifyVector = T)
annot <- fread(js_f$annotation)

whole <- get_seurat(js_f)
setwd(args$out_dir)

## Number of cells before

cells.before <- sapply(whole, function(x)
  dim(GetAssayData(object = x, slot = "counts"))[2])

## NORMALIZATION

whole <- sapply(whole, function(x) SCTransform(
  x,
  ncells=min(100000, ncol(x)),
  vars.to.regress = c("percent.mito"),
  method = "glmGamPoi",
  verbose = T,
  conserve.memory = T
))


## CELL CYCLE SCORING

s.genes <- stringr::str_to_title(cc.genes$s.genes)
g2m.genes <- stringr::str_to_title(cc.genes$g2m.genes)

whole <- sapply(whole, function(x) CellCycleScoring(x, s.features = s.genes,
                                                    g2m.features = g2m.genes,
                                                    set.ident = F,
                                                    assay = 'SCT'))

## INTEGRATION

whole.features <- SelectIntegrationFeatures(object.list = whole, nfeatures = 2000)

whole <- lapply(X = whole, FUN = function(x) {
  x <- RunPCA(x, features = whole.features)
})

whole <- PrepSCTIntegration(object.list = whole, anchor.features = whole.features,
                            verbose = FALSE)
whole.anchors <- FindIntegrationAnchors(object.list = whole, normalization.method = "SCT",
                                        anchor.features = whole.features, verbose = FALSE, reduction = 'rpca')
whole.integrated <- IntegrateData(anchorset = whole.anchors, normalization.method = "SCT",
                                  verbose = FALSE)

## PCA
gc()

whole.integrated <- RunPCA(whole.integrated, verbose = FALSE)


## TSNE

whole.integrated <- RunTSNE(whole.integrated, dims = 1:20, tsne.method = "FIt-SNE",
                            nthreads = 4, max_iter = 2000)
## UMAP

whole.integrated <- RunUMAP(whole.integrated, dims = 1:20)


## CLUSTERING

whole.integrated <- FindNeighbors(object = whole.integrated, dims = 1:20)
whole.integrated <- FindClusters(object = whole.integrated, resolution = c(0.2, 0.4, 0.6, 0.8, 1))

## SAVING: DATASET

write.table(whole.integrated@meta.data, file = 'meta_of_remained_cells.csv', quote = F)

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

idents <- paste0('integrated_snn_res.', c(0.2, 0.4, 0.6, 0.8, 1))
sapply(idents, function(ident) analyze_object(object = whole.integrated, ident = ident))

## Number of cells after

cells.after <- sapply(whole, function(x) length(colnames(x = x)))
cells.diff <- cells.before-cells.after
rbind(cells.before, cells.after, cells.diff)

## Conversion

markers_files <- list.files(pattern = 'markers.tsv', recursive = T, full.names = T)
markers <- lapply(markers_files, function(x) fread(x) %>% dplyr::rename(avg_logFC = 'avg_log2FC'))
names(markers) <- gsub('/|markers', '', stringr::str_extract(markers_files, 'markers/*.*/'))

whole.integrated@meta.data <- whole.integrated@meta.data %>%
  mutate_if(is.character, as.factor)

migrateSeuratObject(whole.integrated,
                    assay="SCT",
                    species='mm',
                    outdir = args$sample_id,
                    public = F,
                    curated = T,
                    markers = markers,
                    generateMarkers = F,
                    generateGMTS = F,
                    name=args$sample_id,
                    token=args$sample_id)
