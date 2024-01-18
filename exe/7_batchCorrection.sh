#!/bin/bash
#SBATCH --job-name=batch
#SBATCH -p veryhimem
#SBATCH --mem=128G       # total memory need
#SBATCH -c 16 #number of cores

loc="/cluster/projects/hansengroup/linyang/PRAD_multomics/seurat"
src="../src"

##-r harmony -t Liger -n the column.names used for batch correction.

Rscript $src/batch_correction.R -i $loc/af_dblt_rm.RDS -o $loc -n "batch" -r true -t true
