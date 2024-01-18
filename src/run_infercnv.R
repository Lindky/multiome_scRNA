###run inferCNV
suppressMessages(library(Seurat))
suppressMessages(library(infercnv))
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input fille path", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="./",
              help="output file path and prefix", metavar="character"),
  make_option(c("-f", "--figure"), type="logical", default="FALSE",
              help="if generate the infercnv plots", metavar="character"),
  make_option(c("-s", "--scale"), type="logical", default="FALSE",
              help="if scaling the data", metavar="character"),
  make_option(c("-t", "--threads"), type="integer", default="16",
              help="number of threads", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

outdir <- opt$outpath

plot <- !opt$figure
print(paste0("infer CNV from sample:", opt$input)) 

infercnv_obj <- readRDS(opt$input)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=outdir,  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             scale_data=opt$scale,
                             analysis_mode="subclusters",
                             output_format="pdf",
                             HMM_type="i6",
                             no_plot=plot,
                             HMM=T,
                             num_threads=opt$threads
)
 
