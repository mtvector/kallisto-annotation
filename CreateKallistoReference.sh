#!/bin/bash

#conda create -n kallisto
#conda activate kallisto
#conda install -c bioconda bustools
#conda install -c bioconda kallisto
#conda install -c bioconda pybedtools
#pip install pandas 

source ~/.bashrc
echo $PATH
which python3
#conda list
#pip list
export PATH="/wynton/home/ye/mschmitz1/utils/miniconda3/bin:$PATH"
conda activate kallisto
which python3
SCRIPTPATH=$1
REFPATH=$2
FA=$3
GTF=$4

cd $REFPATH

#Generating cdna and intron for kallisto index
echo "gffread -w cDNA.fa -g ${FA}.fa ${GTF}.gtf"
gffread -w cDNA.fa -g ${FA}.fa ${GTF}.gtf
echo "python3 ${SCRIPTPATH}/GTF2BED.py ${GTF}.gtf"
python3 ${SCRIPTPATH}/GTF2BED.py ${GTF}.gtf

cut -f1,2 ${FA}.fa.fai | sort -k1,1 -k2,2n  > chromSizes.txt
awk 'OFS="\t" {print $1, "0", $2}' chromSizes.txt | sort -k1,1 -k2,2n > chromSizes.bed
sort -k1,1 -k2,2n ${GTF}.bed > tx.sorted.bed
head -1 tx.sorted.bed

bedtools complement -i tx.sorted.bed -g chromSizes.txt > intergenic_sorted_intermediate.bed
awk 'OFS="\t"{print $1, $2, $3, "--"}' intergenic_sorted_intermediate.bed > intergenic_sorted.bed
head -1 intergenic_sorted.bed

cat ${GTF}.gtf | awk '$3 == "exon" {print $0}' | sort -k1,1 -k4,4n -k5,5n > exon_sorted.gtf
python3 $SCRIPTPATH/GTF2BED.py exon_sorted.gtf
awk 'OFS="\t"{print $1, $2, $3,$4}' exon_sorted.bed > exon_sorted.final.bed

bedtools complement -i <(cat exon_sorted.final.bed intergenic_sorted.bed | sort -k1,1 -k2,2n) -g chromSizes.txt > intron_sorted.bed
cat ${GTF}.gtf | awk '$3 == "transcript" {print $0}' | sort -k1,1 -k4,4n -k5,5n > transcript_sorted.gtf

python3 $SCRIPTPATH/GTF2BED.py transcript_sorted.gtf
#Need bedtools for finding overlap for introns, conflict in kallisto env
python3 $SCRIPTPATH/IntronTranscriptBedOverlap.py intron_sorted.bed transcript_sorted.bed
#conda activate kallisto

cat ${GTF}.gtf | ~/utils/t2g.py > tr2g.txt

bedtools getfasta -name -fo introns.fa -fi ${FA}.fa -bed intron_named.bed
head -2 introns.fa

#Fix the headers on the fasta files
python3 $SCRIPTPATH/CorrectFaHeader.py ${GTF}.gtf cDNA.fa cDNA.correct_header.fa  
python3 $SCRIPTPATH/CorrectFaHeader.py ${GTF}.gtf introns.fa introns.correct_header.fa -I                      

head -1 introns_t2g.txt
head -1 cDNA_t2g.txt

cat introns.correct_header.fa | awk '/^>/ {print $0}' | tr ":" " " | tr -d ">" | awk '{print $1}' > introns_transcripts.txt
cat introns_transcripts.txt | awk '{print $0"."NR"-I"}' > introns_transcripts.to_capture.txt

cat cDNA.correct_header.fa | awk '/^>/ {print $0}'| tr -d ">" | awk '{print $1}' > cDNA_transcripts.txt
cat cDNA_transcripts.txt | awk '{print $0"."NR}' > cDNA_transcripts.to_capture.txt
head -1 cDNA_transcripts.to_capture.txt

cat cDNA.correct_header.fa introns.correct_header.fa > cDNA_introns.fa
cat cDNA_t2g.txt introns_t2g.txt > cDNA_introns_t2g.txt

#Make files for intron exon together capture (better for cellbender)
cat cDNA_introns_t2g.txt | awk '{if ($1 ~ /-I/) {print $1, $2"__I", $3"__I"} else {print $1, $2, $3}}' > cDNA_introns_t2g.markintron.txt
cat cDNA_transcripts.txt > cDNA_introns_transcripts.txt
cat introns_transcripts.txt >> cDNA_introns_transcripts.txt

head -1 cDNA_introns.fa
head -1 cDNA_introns_t2g.txt

kallisto index -i cDNA.idx -k 31 cDNA.fa
kallisto index -i cDNA_introns.idx -k 31 cDNA_introns.fa
