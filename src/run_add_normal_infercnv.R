###run inferCNV

suppressMessages(library(Seurat))
suppressMessages(library(infercnv))
suppressMessages(library(optparse))
suppressMessages(library(stringr))
suppressMessages(library(stringi))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input fille path", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="./",
              help="output file path and prefix", metavar="character"),
  make_option(c("-a", "--anno"), type="character", default="./",
              help="cell annotation file", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default="./",
              help="sampleID", metavar="character"),
  make_option(c("-l", "--label"), type="character", default="./",
              help="the name of annotation column", metavar="character"),
  make_option(c("-m", "--norm"), type="character", default="./",
              help="public", metavar="character"),
  make_option(c("-n", "--scale"), type="logical", default="./",
              help="whether apply normalization", metavar="character"),
  make_option(c("-c", "--control"), type="character", default="./",
              help="label for reference cells", metavar="character"),
  make_option(c("-r", "--ref"), type="character", default="./",
              help="a file indicates the gene location on the chromosome", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

indir <- opt$input
outdir <- opt$outpath
gene_order <- opt$ref
anno_file <- opt$anno
control_label <- opt$control

control_label <- str_split_fixed(control_label, ",", 4)
control_label <- as.character(control_label)
control_label <- as.character(stri_remove_empty(control_label, na_empty = FALSE))

print("reference cells: ")
print(control_label)

data <- readRDS(indir)

sam <- paste0(opt$sample, "_rna")

##extracting the subset samples including the control cells and query cells
tmp_meta <- data@meta.data
#control_cellid <- rownames(tmp_meta[tmp_meta[[opt$label]] %in% control_label,])
sample_id <- rownames(tmp_meta[tmp_meta[["orig.ident"]] == sam,])

#data_id <- c(sample_id, control_cellid)
#data_id <- data_id[!duplicated(data_id)]

data <- data[, which(colnames(data) %in% sample_id)]

print(paste0("query sample sizes:", length(sample_id)))

count_matrix <- data@assays$RNA@counts[, colnames(data)]
#print(paste0("control sample sizes:", length(control_cellid)))
#outdir <- "/Users/linyang/Documents/Hansen_lab/multi-omics/seurat/combine/inferCNV"
norm <- readRDS(opt$norm)
tmp_count <- as.data.frame(count_matrix)
tmp_count <- merge(tmp_count, norm, by.x = 0, by.y = 1, all = FALSE)
rownames(tmp_count) <- tmp_count[[1]]
tmp_count <- tmp_count[-1]
count_matrix <- as.sparse(tmp_count)
###maglinant vs normal 
#load the cell annotation file, use the immune cell as normal control 
ann <- read.delim(anno_file, header = FALSE)

ann <- ann[ann[[1]] %in% colnames(count_matrix),]
print(paste0("subset samples sizes:", nrow(ann)))

bsname <- basename(anno_file)
tmp_ann <- sub(bsname, "tmp_label.txt", anno_file)
write.table(ann, file = tmp_ann, row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

####load the gene order reference 
#gene <- read.delim(gene_order)
print(paste0("infer CNV from sample:", opt$sample)) 
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=count_matrix,
                                    annotations_file=tmp_ann,
                                    delim="\t",
                                    gene_order_file=gene_order,
                                    ref_group_names=c(control_label)) 


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=outdir,  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             scale_data=opt$scale,
                             analysis_mode="subclusters",
                             output_format="pdf",
                             HMM_type="i6",
                             HMM=T,
                             num_threads = 16
)
 
