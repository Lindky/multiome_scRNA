suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input fille path", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="./",
              help="output file path and prefix ", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default="./",
              help="sample id", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

indir <- opt$input
outdir <- opt$outpath
sam <- opt$sample

asam <- paste0(sam,"_atac")
rsam <- paste0(sam,"_rna")

data <- Read10X(data.dir = indir)
pbmc <- CreateSeuratObject(counts = data$`Gene Expression`, project = rsam, min.cells = 3, min.features = 200)
#pbmc.atac <- CreateSeuratObject(counts = data$Peaks, assay = "ATAC", project = asam)

saveRDS(pbmc, file = paste0(outdir, "/RNA.RDS"))
#saveRDS(pbmc.atac, file = paste0(outdir, "/atac.RDS"))



