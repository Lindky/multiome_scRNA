suppressMessages(library(infercnv))
suppressMessages(library(optparse))
suppressMessages(library(stringr))
suppressMessages(library(stringi))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input fille path", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="./",
              help="output file path and prefix", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

data <- readRDS(opt$input)
plot_cnv(data)
