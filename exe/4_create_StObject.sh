#!/bin/bash
#SBATCH --job-name=ct_seurat
#SBATCH -p himem
#SBATCH --mem=32G       # total memory need
#SBATCH -c 8 #number of cores

inputdir="/cluster/projects/hansengroup/linyang/PRAD_multomics/third_batch"
outdir="/cluster/projects/hansengroup/linyang/PRAD_multomics"

src="../src"

mkdir $outdir/seurat/

for sam in $(cat ${inputdir}/samlist.txt)
do
  echo "${sam}"
  mkdir $outdir/seurat/${sam}
  
  Rscript $src/create_srObject.R -i $inputdir/$sam/outs/filtered_feature_bc_matrix \
          -o $outdir/seurat/${sam} -s $sam

done
