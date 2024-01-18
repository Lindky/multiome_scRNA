#!/bin/bash
#SBATCH --job-name=doublet
#SBATCH -p veryhimem
#SBATCH --mem=128G       # total memory need
#SBATCH -c 16 #number of cores

loc="/cluster/projects/hansengroup/linyang/PRAD_multomics/seurat"
src="../src"
Rscript $src/doubletfinder.R -i $loc/bf_dblt_data.RDS -o $loc -s $loc/summary_tab.RDS \
        -c $loc/batch2_3_clinical_info.csv
