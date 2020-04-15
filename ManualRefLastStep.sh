#!/bin/bash
#$ -o ~/log
#$ -e ~/log
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l mem_free=48G
#$ -l h_rt=6:00:00

source ~/.bashrc
conda activate kallisto
REFPATH=/wynton/group/ye/mtschmitz/refdata2/hg38/gencodev33/
FA=hg38
GTF=gencode.v33.annotation

cd $REFPATH
kallisto index -i cDNA.idx -k 31 cDNA.fa
kallisto index -i cDNA_introns.idx -k 31 cDNA_introns.fa


