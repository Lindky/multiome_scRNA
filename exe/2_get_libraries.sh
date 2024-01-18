#!/bin/bash
path="./"
sample_file="samlist.txt"
rna_path="/cluster/projects/hansengroup/linyang/prostate_batch3/20231025_LH00244_0022_A22CYHJLT3_He_Catherine"
atac_path="/cluster/projects/hansengroup/linyang/prostate_batch3/20231101_LH00244_0025_B22CYHVLT3_He_Catherine"

src="../src"

mkdir $path/libraries
for sample in $(cat $path/${sample_file})
do
#in=`expr ${n} + 1`
python $src/make_libraries.py -s ${sample} -o $path/libraries/${sample}.csv -r $rna_path -a $atac_path

done
