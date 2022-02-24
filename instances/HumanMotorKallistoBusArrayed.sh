#$ -V
#$ -S /bin/bash
#$ -e ~/log
#$ -o ~/log
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l mem_free=96G
#$ -l h_rt=200:00:00
#$ -l scratch=200G
#$ -t 1-69:1
#Figure out above number with
#ls *.fastq.gz | awk -F'_S[0-9]+' '{print $1}' | uniq | wc -l

JOB_DIR=/wynton/scratch/mtschmitz/humanfastqpool
source "/wynton/home/ye/mschmitz1/.bashrc"
conda activate kallisto
mkdir -p $JOB_DIR
SGE_TASK_ID=`expr $SGE_TASK_ID - 1`
#gunzip ~/utils/cellranger-3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/*
threads=1
genome=/wynton/group/ye/mtschmitz/refdata2/hg38
echo $SGE_TASK_ID
cd $JOB_DIR
files=($(ls $JOB_DIR/*.fastq.gz | awk -F'_S[0-9]+' '{print $1}' | uniq))
file="${files[$SGE_TASK_ID]}"
echo $file
#samplename=$(echo $filebase | awk -F "_S[0-9]" '{print $NR}')
#samplefiles=$(ls -1 $JOB_DIR/${samplename}*.fastq* | xargs -n 1 basename |  paste -sd " " -)

echo ${file}_kOut/output.correct.sort.bus

if [ -r ${file}_kOut/output.correct.sort.bus ]
then
echo "Done"
else
readlen=$(zcat ${file}*R1* | head -4 | awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}')
echo $readlen
if [[ $readlen == 26 ]];
then
echo "kallisto bus -t $threads -i ${genome}/cDNA_introns.idx -o ${file}_kOut/ -x 10xv2 ${file}*.fastq*"
kallisto bus -t $threads -i ${genome}/cDNA_introns.idx -o ${file}_kOut/ -x 10xv2 ${file}*.fastq*
cd ${file}_kOut/
mkdir -p cDNA_capture/ introns_capture/ spliced/ unspliced/ tmp/
bustools correct -w ~/utils/cellranger-3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt -p output.bus | bustools sort -T tmp/ -o output.correct.sort.bus -t $threads -
else
echo "kallisto bus -t $threads -i ${genome}/cDNA_introns.idx -o ${file}_kOut/ -x 10xv3 ${file}*.fastq*"
kallisto bus -t $threads -i ${genome}/cDNA_introns.idx -o ${file}_kOut/ -x 10xv3 ${file}*.fastq*
cd ${file}_kOut/
mkdir -p cDNA_capture/ introns_capture/ spliced/ unspliced/ tmp/
bustools correct -w ~/utils/cellranger-3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt -p output.bus | bustools sort -T tmp/ -o output.correct.sort.bus -t $threads -
fi
fi
cd  ${file}_kOut/
bustools capture -o cDNA_capture/ -c $genome/cDNA_transcripts.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus
bustools capture -o introns_capture/ -c $genome/introns_transcripts.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus
bustools count -o unspliced/u -g $genome/cDNA_introns_t2g.txt -e cDNA_capture/split.ec -t transcripts.txt --genecounts cDNA_capture/split.bus
bustools count -o spliced/s -g $genome/cDNA_introns_t2g.txt -e introns_capture/split.ec -t transcripts.txt --genecounts introns_capture/split.bus
cd $JOB_DIR

#mkdir -p /wynton/group/ye/mtschmitz/kallistoouts
#cp -r $JOB_DIR/*kOut /wynton/group/ye/mtschmitz/kallistoouts
