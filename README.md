# Overview
This is a pipeline to process single-cell RNA data from [10X multiome](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression). 

**Note:** If your data is single-cell RNA **only** from the 10X, you have to run cellranger by yourself.  
After running cellranger, you can start from **4_create_StObject.sh**

# Tutorial
## Pre-requisites:
In order to run this pipeline, you have to install the following packages:

**Seurat**: https://github.com/satijalab/seurat  
**optparse**: https://www.rdocumentation.org/packages/optparse/versions/1.7.3  
**SeuratWrappers**: https://github.com/satijalab/seurat-wrappers  
**DoubletFinder**: https://github.com/chris-mcginnis-ucsf/DoubletFinder  


(optional): if you want to do batch correction using harmony or Liger   
**harmony**: https://github.com/immunogenomics/harmony  
**rliger**: https://github.com/welch-lab/liger


## Run the pipeline 

The the exection file is the file with the ".sh" extension in **./exe** folder.  
The number of each file indicating the order of analysis.


The scripts used in the exection file are saved under the **./src** folder.



