suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(DropletUtils))
suppressMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument('--rda',
                    type = "character",
                    help = 'Path to directory with barcodes / features / matrix files')
parser$add_argument('--out_dir',
                    type = "character",
                    help = 'Path to output directory')
parser$add_argument('--json',
                    type = "character",
                    help = 'Path to json file with target clusters')
parser$add_argument('--name',
                    type = "character",
                    help = 'Token corresponds to json key')
parser$add_argument('--annot',
                    type = "character",
                    help = 'path to metadata per sample')
args <- parser$parse_args()
print(args)

get_data <- function(path) {
  load(path)
  whole.integrated
}

filter_data <- function(object, js_f, name, annot) {
  if (!is.null(js_f[[name]])) {
    cells_to_remove <- object@meta.data %>%
      filter(!!as.symbol(names(js_f[[name]])) %in%
               strsplit(js_f[[name]][1][[1]], " ")[[1]]) %>%
      rownames()
    object <- subset(object, cells = cells_to_remove, invert = T)
  }
  annot <- fread(annot)
  object@meta.data <- object@meta.data %>%
    dplyr::select(grep('RNA|mito|sample|Feature|Count',
                        colnames(object@meta.data), value = T))
  if (name == 'GSE163005') {
    cells_to_remove <- object@meta.data %>% 
      filter(!grepl('MS19270|MS49131|MS58637|MS71658', sample)) %>% 
      rownames()
    object <- subset(object, cells = cells_to_remove)
  }
  object$study <- annot$study[match(object$sample, annot$sample)]
  object$sex <- annot$sex[match(object$sample, annot$sample)]
  object$organ <- annot$organ[match(object$sample, annot$sample)]
  object$disease <- annot$disease[match(object$sample, annot$sample)]
  object
}

write_data <- function(object, out_dir) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }
  write.table(object@meta.data, file=sprintf('%s/meta_data.csv', out_dir),
              sep=',', row.names = T, quote = F)
  write10xCounts(sprintf('%s/data.h5', out_dir), object@assays$RNA@counts,
                 gene.symbol=rownames(object@assays$RNA@counts),
                 barcodes=colnames(object@assays$RNA@counts))
}

js_f <- jsonlite::read_json(args$json, simplifyVector = T)

print(js_f)

object <- get_data(args$rda)

object <- filter_data(object, js_f, args$name, args$annot)

write_data(object, args$out_dir)

