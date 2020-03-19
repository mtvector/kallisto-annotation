#!/bin/bash
#$ -o ~/log
#$ -e ~/log
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l mem_free=4G
#$ -l h_rt=6:00:00

source ~/.bashrc
conda activate kallisto

REFPATH=/wynton/group/ye/mtschmitz/refdata2/copyrhemac10/
FA=rheMac10
GTF=rheMac10_cat_hg38_V32_Lv3_v2.corrected

echo "$PWD $REFPATH $FA $GTF"
$PWD/CreateKallistoReference.sh $PWD $REFPATH $FA $GTF
