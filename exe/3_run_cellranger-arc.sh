#!/bin/bash
#SBATCH --job-name=cr-arc
#SBATCH -p veryhimem
#SBATCH --mem=256G       # total memory need
#SBATCH -c 32 #number of cores

#SBATCH --array=0-5%5


module load cellranger-arc/2.0.2
module load igenome-human/hg38

INPUT=( $(cat ./samlist.txt) )
ref="/cluster/tools/data/commondata/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
library_loc="./libraries"

func () {
  sample=$1 
  echo "runing the sample: ${sample}..."
  
  libraryFile=$library_loc/${sample}.csv

  cellranger-arc count --id="${sample}" \
                       --reference=$ref \
                       --libraries=$libraryFile \
                       --localcores=32 \
                       --localmem=256
}

func "${INPUT[$SLURM_ARRAY_TASK_ID]}"
