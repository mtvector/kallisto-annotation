#!/bin/bash

#conda create -n kallisto
#conda activate kallisto
#conda install -c bioconda bustools
#conda install -c bioconda kallisto
#conda install -c bioconda pybedtools
#pip install pandas 

source activate kallisto

SCRIPTPATH=$1
REFPATH=$2
FA=$3
GTF=$4

cd $REFPATH

#Generating cdna and intron for kallisto index
echo "gffread -w ${FA}.cDNA.fa -g  ${FA}.fa ${GTF}.gtf"
gffread -w cDNA.fa -g  ${FA}.fa ${GTF}.gtf
python ${SCRIPTPATH}/GTF2BED.py ${GTF}.gtf

cut -f1,2 ${FA}.fa.fai | sort -k1,1 -k2,2n  > chromSizes.txt
awk 'OFS="\t" {print $1, "0", $2}' chromSizes.txt | sort -k1,1 -k2,2n > chromSizes.bed
sort -k1,1 -k2,2n ${GTF}.bed > tx.sorted.bed
head -1 tx.sorted.bed

bedtools complement -i tx.sorted.bed -g chromSizes.txt > intergenic_sorted_intermediate.bed
awk 'OFS="\t"{print $1, $2, $3, "--"}' intergenic_sorted_intermediate.bed > intergenic_sorted.bed
head -1 intergenic_sorted.bed

cat ${GTF}.gtf | awk '$3 == "exon" {print $0}' | sort -k1,1 -k4,4n -k5,5n > exon_sorted.gtf
python $PATH/GTF2BED.py/GTF2BED.py exon_sorted.gtf
awk 'OFS="\t"{print $1, $2, $3,$4}' exon_sorted.bed > exon_sorted.final.bed

bedtools complement -i <(cat exon_sorted.final.bed intergenic_sorted.bed | sort -k1,1 -k2,2n) -g chromSizes.txt > intron_sorted.bed
cat ${GTF}.gtf | awk '$3 == "transcript" {print $0}' | sort -k1,1 -k4,4n -k5,5n > transcript_sorted.gtf

python $PATH/GTF2BED.py/GTF2BED.py transcript_sorted.gtf
#Need bedtools for finding overlap for introns, conflict in kallisto env
#conda activate biopython
python $PATH/GTF2BED.py/IntronTranscriptBedOverlap.py intron_sorted.bed transcript_sorted.bed
#conda activate kallisto

cat ${GTF}.gtf | ~/utils/t2g.py > tr2g.txt

bedtools getfasta -name -fo introns.fa -fi ${FA}.fa -bed intron_named.bed
head -2 introns.fa

#Fix the headers on the fasta files
python $PATH/GTF2BED.py/CorrectFaHeader.py ${GTF}.gtf cDNA.fa cDNA.correct_header.fa  
python $PATH/GTF2BED.py/CorrectFaHeader.py ${GTF}.gtf introns.fa introns.correct_header.fa -I                      

head -1 introns_t2g.txt
head -1 cDNA_t2g.txt

cat introns.correct_header.fa | awk '/^>/ {print $0}' | tr ":" " " | tr -d ">" | awk '{print $1}' > introns_transcripts.txt
cat introns_transcripts.txt | awk '{print $0"."NR"-I"}' > introns_transcripts.to_capture.txt
#awk 'NR==FNR{a[$1]=$2; b[$1]=$3;next} {$2=a[$1];$3=b[$1]} 1' tr2g.txt introns_transcripts.txt > introns_t2g.txt

#awk '{print ">"$1"."NR"-I"" gene_id:"$2" gene_name:"$3}' introns_t2g.txt > introns_fasta_header.txt
#awk -v var=1 'FNR==NR{a[NR]=$0;next}{ if ($0~/^>/) {print a[var], var++} else {print $0}}' introns_fasta_header.txt introns.fa > introns.correct_header.fa

cat cDNA.correct_header.fa | awk '/^>/ {print $0}'| tr -d ">" | awk '{print $1}' > cDNA_transcripts.txt
cat cDNA_transcripts.txt | awk '{print $0"."NR}' > cDNA_transcripts.to_capture.txt
head -1 cDNA_transcripts.to_capture.txt

#awk 'NR==FNR{a[$1]=$2; b[$1]=$3;next} {$2=a[$1];$3=b[$1]} 1' tr2g.txt cDNA_transcripts.txt > cDNA_t2g.txt
#awk '{print ">"$1"."NR" gene_id:"$2" gene_name:"$3}' cDNA_t2g.txt > cDNA_fasta_header.txt
#awk -v var=1 'FNR==NR{a[NR]=$0;next}{ if ($0~/^>/) {print a[var], var++} else {print $0}}' cDNA_fasta_header.txt $cDNA_fa > cDNA.correct_header.fa
#head -1 cDNA.correct_header.fa

cat cDNA.correct_header.fa introns.correct_header.fa > cDNA_introns.fa
cat cDNA_t2g.txt introns_t2g.txt > cDNA_introns_t2g.txt
head -1 cDNA_introns.fa
head -1 cDNA_introns_t2g.txt

kallisto index -i cDNA.idx -k 31 cDNA.fa
kallisto index -i cDNA_introns.idx -k 31 cDNA_introns.fa
