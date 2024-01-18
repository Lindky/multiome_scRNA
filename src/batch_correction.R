suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(optparse))
suppressMessages(library(harmony))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(rliger))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input fille path", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="./",
              help="output file path and prefix ", metavar="character"),
  make_option(c("-r", "--harmony"), type="logical", default="False",
              help="apply harmony", metavar="character"),
  make_option(c("-n", "--nor"), type="character", default="None",
              help="batch for batch correction", metavar="character"),
  make_option(c("-t", "--liger"), type="logical", default="False",
              help="apply liger", metavar="character")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

indir <- opt$input
outdir <- opt$outpath

print("starting...")
set.seed(100)

mgh.combine <- readRDS(indir)

##cluster without doing batch correction 
mgh.combine <- SCTransform(mgh.combine, verbose = TRUE)
mgh.combine <- RunPCA(mgh.combine, features = VariableFeatures(object = mgh.combine,
                                                               assay = "SCT"), assay = "SCT")
mgh.combine <- FindNeighbors(mgh.combine, dims = 1:30, verbose = TRUE, assay = "SCT",
                     reduction = "pca", graph.name = c("NN.SCT", "SNN.SCT"))
mgh.combine <- FindClusters(mgh.combine, verbose = TRUE, graph.name = "SNN.SCT")
mgh.combine <- RunUMAP(mgh.combine, reduction = "pca", dims = 1:30, verbose = TRUE,
               reduction.name = "UMAP.SCT", assay = "SCT")

##convert the batch colum as factors
meta <- mgh.combine@meta.data
meta[[opt$nor]] <- paste0(meta[[opt$nor]], "_batch")
meta[[opt$nor]] <- as.factor(meta[[opt$nor]])
mgh.combine@meta.data <- meta


if(opt$harmony){

   print("Batch correction using harmony...")
   print(paste0("Using covariant: ", opt$nor))

   mgh.combine <- mgh.combine %>% RunHarmony(opt$nor)
 
   mgh.combine <- RunUMAP(mgh.combine, reduction = "harmony", reduction.name = "UMAP.harmony",dims = 1:30)
   mgh.combine <- FindNeighbors(mgh.combine, reduction = "harmony", dims = 1:30,
                                graph.name = c("NN.harmony", "SNN.harmony"))
   mgh.combine <- FindClusters(mgh.combine, verbose = TRUE, graph.name = "SNN.harmony")

   print("Done")
}

ifnb <- mgh.combine

if(opt$liger){
  print("Batch correction using iNMF...")
  print(paste0("Using covariant: ", opt$nor))

  DefaultAssay(object = ifnb) <- "RNA"
  ifnb <- NormalizeData(ifnb)
  ifnb <- FindVariableFeatures(ifnb)
  ifnb <- ScaleData(ifnb, split.by = opt$nor, do.center = FALSE)
  
  ifnb <- RunOptimizeALS(ifnb, k = 30, lambda = 5, split.by = opt$nor)
  ifnb <- RunQuantileNorm(ifnb, split.by = opt$nor)

  ifnb <- FindNeighbors(ifnb, reduction = "iNMF", dims = 1:30,
                       graph.name = c("NN.inmf", "SNN.inmf"))
  ifnb <- FindClusters(ifnb, verbose = TRUE, graph.name = "SNN.inmf")
  ifnb <- RunUMAP(ifnb, dims = 1:ncol(ifnb[["iNMF"]]), reduction = "iNMF", reduction.name = "UMAP.iNMF")

  print("Done")
}

saveRDS(ifnb, file = paste0(outdir, "/batch_correction.RDS"))







