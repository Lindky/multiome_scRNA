suppressMessages(library(Seurat))
suppressMessages(library(optparse))
suppressMessages(library(DoubletFinder))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input fille path", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="./",
              help="output file path and prefix ", metavar="character"),
  make_option(c("-s", "--summary"), type="character", default="False",
              help="summary table", metavar="character"),
  make_option(c("-c", "--clinical"), type="character", default="False",
              help="clinical table", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

indir <- opt$input
outdir <- opt$outpath

print("starting...")
set.seed(100)

summary_table <- readRDS(opt$summary)
mgh_post <- readRDS(indir)

print("Identify doublets using DoubletFinder ...")

mgh_doubletfinder_total <- lapply(mgh_post, function(x){
 
  sam <- head(x@meta.data$orig.ident,1)
  print(" ")
  print(" ")
  print(paste0("processing the sample: ", sam))

  sweep <- paramSweep_v3(x, PCs = 1:30, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  find_pk <- bcmvn$pK[bcmvn$BCmetric==max(bcmvn$BCmetric)]
  # levels(find_pk)
  find_pk_num <- levels(find_pk)[as.numeric(find_pk)]
  find_pk_num <- as.numeric(find_pk_num)
  
  nExp_poi <- round(0.03*nrow(x@meta.data))
  
  #Homotypic Doublet Proportion Estimate
  annotations <- x@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  biopsy_doublet <- doubletFinder_v3(x, PCs = 1:30, pN = 0.25, pK = find_pk_num, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  
  reuse_col <- colnames(biopsy_doublet@meta.data)[ncol(biopsy_doublet@meta.data)-1]
  
  biopsy_doublet_heterotypic <- doubletFinder_v3(biopsy_doublet, PCs = 1:30, pN = 0.25, pK = find_pk_num, nExp = nExp_poi.adj, reuse.pANN = reuse_col, sct = TRUE)
  
  print("Done...") 
  return(biopsy_doublet_heterotypic)
})

#saveRDS(mgh_doubletfinder_total, file = paste0(outdir, "/check_doublet_data.RDS"))

summary_dblt <- NULL
for (sam in summary_table$sample) {
  tmp.data <- mgh_doubletfinder_total[[sam]]
  tmp.meta <- tmp.data@meta.data
  
  colnames(tmp.meta)[ncol(tmp.meta)] <- "doublet"
  n = ncol(tmp.meta) -1 
  tmp.meta <- tmp.meta[-n]
  
  tmp.data@meta.data <- tmp.meta
     
  doublet_plot <- DimPlot(tmp.data, group.by = "doublet")

  pdf(paste0(outdir,"/",sam, "/doublet_umap.pdf"))
  print(doublet_plot)
  dev.off()

  dblt_ratio <- nrow(tmp.meta[tmp.meta$doublet == "Doublet",])/nrow(tmp.meta)
  tab <- data.frame(sample = sam, dblt_ratio = dblt_ratio)
  summary_dblt <- rbind(summary_dblt, tab)
  mgh_doubletfinder_total[[sam]] <- tmp.data
}

summary_table <- merge(summary_table, summary_dblt, by = "sample")
write.csv(summary_table, file = paste0(outdir, "/qc_summary.csv"), quote = F, row.names = F)

print("start merging samples...")

mgh.combine <- merge(mgh_doubletfinder_total[[1]], mgh_doubletfinder_total[[2]])

num <- length(mgh_doubletfinder_total)
if(num >= 3) {
  for (n in 3:num) {
    mgh.combine <- merge(mgh.combine, mgh_doubletfinder_total[[n]])
  }
}

meta <- mgh.combine@meta.data
remove_id <- grep("pANN_", colnames(meta))
meta <- meta[-remove_id]
mgh.combine@meta.data <- meta

clinical_tab <- read.csv(opt$clinical)
clinical_tab[[1]] <- paste0(clinical_tab[[1]], "_rna")

meta$num <- 1:nrow(meta)
meta <- merge(meta, clinical_tab, by.x = "orig.ident", by.y = 1, all = TRUE)
meta <- meta[order(meta$num, decreasing = F),]

meta_terms <- colnames(clinical_tab)[-1]
for (m in meta_terms) {
  mgh.combine[[m]] <- meta[[m]]
}
saveRDS(mgh.combine, file = paste0(outdir, "/merged_data.RDS"))

print("Done")

