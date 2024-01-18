#!/bin/bash

data_loc="/cluster/projects/hansengroup/linyang/prostate_batch3/20231101_LH00244_0025_B22CYHVLT3_He_Catherine"
path="./"
#generate the id text file first 
ls $data_loc > ${path}/id.txt

for name in $(cat ./id.txt)
do
echo "${name%%_S*}"

done >> tmp.txt  

sort -u tmp.txt > samlist.txt
rm tmp.txt
rm id.txt

