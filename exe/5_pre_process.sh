#!/bin/bash
#SBATCH --job-name=seurat
#SBATCH -p veryhimem
#SBATCH --mem=128G       # total memory need
#SBATCH -c 32 #number of cores

loc="/cluster/projects/hansengroup/linyang/PRAD_multomics/seurat"
src="../src"
Rscript $src/dimentional_rd.R -i $loc -o $loc
