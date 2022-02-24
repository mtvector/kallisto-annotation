#$ -V
#$ -S /bin/bash
#$ -e ~/log
#$ -o ~/log
#$ -cwd
#$ -j y
#$ -pe smp 2
#$ -l mem_free=96G
#$ -l h_rt=200:00:00
#$ -l scratch=200G

JOB_DIR=/wynton/scratch/mtschmitz/humanfastqpool
source "/netapp/home/mschmitz1/.bashrc"
conda activate kallisto
mkdir -p $JOB_DIR

gunzip ~/utils/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/*
threads=2
genome=/wynton/group/ye/mtschmitz/refdata2/hg38

#cp -r /ye/yelabstore2/mtschmitz/seq/refdata/$genome $JOB_DIR
cd $JOB_DIR
files=($JOB_DIR/*.f*)
#for ((i=${#files[@]}-1; i>=0; i--)); do
for ((i=0; i<=${#files[@]}-1; i++)); do

file="${files[$i]}"
filebase=${file##*/}
samplename=$(echo $filebase | awk -F "_S[0-9]" '{print $NR}')
samplefiles=$(ls -1 $JOB_DIR/${samplename}*.fastq* | xargs -n 1 basename |  paste -sd " " -)

echo $filebase
echo $samplefiles

#if [ -r $JOB_DIR/${samplename}_kOut/output.bus ]
if [ -r $JOB_DIR/${samplename}_kOut/output.correct.sort.bus ]
then
echo "Done"
else
readlen=$(zcat $JOB_DIR/${samplename}*R1* | head -4 | awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}')
echo $readlen
if [[ $readlen == 26 ]];
then
kallisto bus -t $threads -i ${genome}/cDNA_introns.idx -o ${samplename}_kOut/ -x 10xv2 $JOB_DIR/${samplename}*.fastq*
cd ${samplename}_kOut/
mkdir -p cDNA_capture/ introns_capture/ spliced/ unspliced/ tmp/
bustools correct -w ~/utils/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/737K-august-2016.txt -p output.bus | bustools sort -T tmp/ -o output.correct.sort.bus -t $threads -
else
kallisto bus -t $threads -i ${genome}/cDNA_introns.idx -o ${samplename}_kOut/ -x 10xv3 $JOB_DIR/${samplename}*.fastq*
cd ${samplename}_kOut/
mkdir -p cDNA_capture/ introns_capture/ spliced/ unspliced/ tmp/
bustools correct -w ~/utils/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/3M-february-2018.txt -p output.bus | bustools sort -T tmp/ -o output.correct.sort.bus -t $threads -
fi
fi
cd  $JOB_DIR/${samplename}_kOut/
bustools capture -o cDNA_capture/split.bus-c $genome/cDNA_transcripts.txt -s -e matrix.ec -t transcripts.txt output.correct.sort.bus
bustools capture -o introns_capture/split.bus -c $genome/introns_transcripts.txt -s -e matrix.ec -t transcripts.txt output.correct.sort.bus
bustools count -o unspliced/u -g $genome/cDNA_introns_t2g.txt -e cDNA_capture/split.ec -t transcripts.txt --genecounts cDNA_capture/split.bus
bustools count -o spliced/s -g $genome/cDNA_introns_t2g.txt -e introns_capture/split.ec -t transcripts.txt --genecounts introns_capture/split.bus
cd $JOB_DIR
done

mkdir -p /wynton/group/ye/mtschmitz/kallistoouts
cp -r $JOB_DIR/*kOut /wynton/group/ye/mtschmitz/kallistoouts
