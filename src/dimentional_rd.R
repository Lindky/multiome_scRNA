#suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(optparse))
suppressMessages(library(DoubletFinder))
#suppressMessages(library(harmony))
#suppressMessages(library(SeuratWrappers))
#suppressMessages(library(rliger))

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
sample <- list.files(indir)
print(paste0("processing the sample: ", sample))

mgh.data <- sapply(sample, function(i){
  print(file.path(indir,i,"RNA.RDS"))
  d10x <- readRDS(file.path(indir,i,"RNA.RDS"))
  d10x
})

#Quality control 
min.cells <- 10
min.features <- 300

percent.mito.thresh <- 10

summary_table <- NULL
for (s in sample) {
  data <- mgh.data[[s]]
  
  bf_ncell <- ncol(data)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  
  qc_plot <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  pdf(paste0(indir,"/",s, "/pre_qc_rna.pdf"))
  print(qc_plot)
  dev.off()
  
  max_rna <- quantile(data$nFeature_RNA, 0.95)
  max_count <- quantile(data$nCount_RNA, 0.95)
  
  data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < max_rna)
  data <- subset(data, subset = nCount_RNA <= max_count)
  data <- data[,-which(data$percent.mt >= percent.mito.thresh)]
  
  qc_plot <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  pdf(paste0(indir,"/",s, "/af_qc_rna.pdf"))
  print(qc_plot)
  dev.off()
  
  cellbarcode <- colnames(data)
  cellbarcode <- paste0(s, "_", cellbarcode)
  data <- RenameCells(data, new.names = cellbarcode)
 
  mgh.data[[s]] <- data
  summary_tab_tmp <- data.frame(bf_qc = bf_ncell, af_qc = ncol(data))
  summary_table <- rbind(summary_table, summary_tab_tmp)
}
summary_table$sample <- sample

print("perform clustering on each sample...")
#normalization and clustering 
mgh_post <- lapply(mgh.data, function(x){
  
  x <- SCTransform(x, verbose = TRUE)
  x <- RunPCA(x, features = VariableFeatures(object = x, assay = "SCT"), assay = "SCT")
  x <- FindNeighbors(x, dims = 1:30, verbose = TRUE, assay = "SCT", reduction = "pca",   graph.name = "SNN.SCT")
  x <- FindClusters(x, verbose = TRUE, graph.name = "SNN.SCT")
  x <- RunUMAP(x, dims = 1:30, verbose = TRUE, reduction.name = "UMAP.SCT", assay = "SCT")
  ##### add SCTransform, PCA, RunTSNE for each sample, do the quality plot for each sample and save them in the folder, and we do the doublet remove
  return(x)
})

saveRDS(summary_table, file = paste0(outdir, "/summary_tab.RDS"))
saveRDS(mgh_post, file = paste0(outdir, "/bf_dblt_data.RDS"))

#print("Identify doublets using DoubletFinder ...")
#
#mgh_doubletfinder_total <- lapply(mgh_post, function(x){
# 
#  sam <- head(x@meta.data$orig.ident,1)
#  print(" ")
#  print(" ")
#  print(paste0("processing the sample: ", sam))
#
#  sweep <- paramSweep_v3(x, PCs = 1:30, sct = TRUE)
#  sweep.stats <- summarizeSweep(sweep, GT = FALSE)
#  bcmvn <- find.pK(sweep.stats)
#  find_pk <- bcmvn$pK[bcmvn$BCmetric==max(bcmvn$BCmetric)]
#  # levels(find_pk)
#  find_pk_num <- levels(find_pk)[as.numeric(find_pk)]
#  find_pk_num <- as.numeric(find_pk_num)
#  
#  nExp_poi <- round(0.03*nrow(x@meta.data))
#  
#  #Homotypic Doublet Proportion Estimate
#  annotations <- x@meta.data$seurat_clusters
#  homotypic.prop <- modelHomotypic(annotations)
#  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#  
#  biopsy_doublet <- doubletFinder_v3(x, PCs = 1:30, pN = 0.25, pK = find_pk_num, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
#  
#  reuse_col <- colnames(biopsy_doublet@meta.data)[ncol(biopsy_doublet@meta.data)-1]
#  
#  biopsy_doublet_heterotypic <- doubletFinder_v3(biopsy_doublet, PCs = 1:30, pN = 0.25, pK = find_pk_num, nExp = nExp_poi.adj, reuse.pANN = reuse_col, sct = TRUE)
#  
#  return(biopsy_doublet_heterotypic)
#})
#
#summary_dblt <- NULL
#for (sam in summary_table$sample) {
#  tmp.data <- mgh_doubletfinder_total[[sam]]
#  tmp.meta <- tmp.data@meta.data
#  
#  colnames(tmp.meta)[ncol(tmp.meta)] <- "doublet"
#  n = ncol(tmp.meta) -1 
#  tmp.meta <- tmp.meta[-n]
#  
#  doublet_plot <- DimPlot(data, group.by = "doublet")
#
#  pdf(paste0(indir,"/",sam, "/doublet_umap.pdf"))
#  print(doublet_plot)
#  dev.off()
#
#  tmp.data@meta.data <- tmp.meta
#  dblt_ratio <- nrow(tmp.meta[tmp.meta$doublet == "Doublet",])/nrow(tmp.meta)
#  tab <- data.frame(sample = sam, dblt_ratio = dblt_ratio)
#  summary_dblt <- rbind(summary_dblt, tab)
#  mgh_doubletfinder_total[[sam]] <- tmp.data
#}
#
#summary_table <- merge(summary_table, summary_dblt, by = "sample")
#write.csv(summary_table, file = paste0(outdir, "/qc_summary.csv"), quote = F, row.names = F)
#
#mgh.combine <- merge(mgh_doubletfinder_total[[1]], mgh_doubletfinder_total[[2]])
#
#num <- length(mgh_doubletfinder_total)
#if(num >= 3) {
#  for (n in 3:num) {
#    mgh.combine <- merge(mgh.combine, mgh_doubletfinder_total[[n]])
#  }
#}
#
#saveRDS(mgh.combine, file = paste0(outdir, "/merged_data.RDS"))

#mgh.combine <- readRDS(indir)
#mgh.combine <- SCTransform(mgh.combine, verbose = TRUE, assay = "Spatial")

#mgh.combine <- RunPCA(mgh.combine, features = VariableFeatures(object = mgh.combine, assay = "SCT"), assay = "SCT")

#Dimensional reduction
#mgh.combine <- FindNeighbors(mgh.combine, dims = 1:30, verbose = TRUE, assay = "SCT", reduction = "pca", graph.name = "SNN.SCT")
#mgh.combine <- FindClusters(mgh.combine, verbose = TRUE, graph.name = "SNN.SCT")
#mgh.combine <- RunUMAP(mgh.combine, dims = 1:30, verbose = TRUE, reduction.name = "UMAP.SCT", assay = "SCT")

#saveRDS(mgh.combine, file = paste0(outdir, "/data.RDS"))


#if(opt$harmony){
#  mgh.combine <- SCTransform(mgh.combine, verbose = TRUE)
#  mgh.combine <- RunPCA(mgh.combine, features = VariableFeatures(object = mgh.combine, assay = "SCT"), assay = "SCT")
#  mgh.combine <- mgh.combine %>% RunHarmony(opt$nor)
 
#  mgh.combine <- RunUMAP(mgh.combine, reduction = "harmony", reduction.name = "UMAP.harmony",dims = 1:30)
#  mgh.combine <- FindNeighbors(mgh.combine, reduction = "harmony", dims = 1:30)
#  mgh.combine <- FindClusters(mgh.combine)
#}

#if(opt$liger){
#  ifnb <- mgh.combine
#  DefaultAssay(object = ifnb) <- "RNA"
#  ifnb <- NormalizeData(ifnb)
#  ifnb <- FindVariableFeatures(ifnb)
#  ifnb <- ScaleData(ifnb, split.by = opt$nor, do.center = FALSE)
  
#  ifnb <- RunOptimizeALS(ifnb, k = 30, lambda = 5, split.by = opt$nor)
#  ifnb <- RunQuantileNorm(ifnb, split.by = opt$nor)

#  ifnb <- FindNeighbors(ifnb, reduction = "iNMF", dims = 1:30)
#  ifnb <- FindClusters(ifnb)
#  ifnb <- RunUMAP(ifnb, dims = 1:ncol(ifnb[["iNMF"]]), reduction = "iNMF", reduction.name = "UMAP.iNMF")
#}

#saveRDS(ifnb, file = paste0(outdir, "/liger.RDS"))


