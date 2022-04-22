#!/bin/bash
#$ -o ~/log
#$ -e ~/log
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l mem_free=300G
#$ -l h_rt=39:00:00

source ~/.bashrc
conda activate kallisto2

REFPATH=/wynton/group/ye/mtschmitz/refdata2/rhemac10/CAT_chang
FA=GCF_003339765.1_Mmul_10_fixed
GTF=Rhesus

echo "$PWD $REFPATH $FA $GTF"
$PWD/CreateKallistoReference.sh $PWD $REFPATH $FA $GTF
