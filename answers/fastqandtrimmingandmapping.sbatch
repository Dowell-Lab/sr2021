#!/bin/bash 
#SBATCH --job-name=fastqcandtrim # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=4     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=01:00:00 # Time limit hrs:min:sec
#SBATCH --output=/Users/maallen3/sread2021/day5/eofiles/fastqcandtrim.%j.out # Standard output
#SBATCH --error=/Users/maallen3/sread2021/day5/eofiles/fastqcandtrim.%j.err # Standard error log


module load fastqc/0.11.5 
module load hisat2/2.1.0
module load samtools/1.8


#The first line in this script tells the computer you want bash (bourne again shell)to run this job.
#the top is the SBATCH information
#the bottom is the command

FILENAMES[0]=chr21Ethan_repA.RNA
FILENAMES[1]=chr21Ethan_repB.RNA
FILENAMES[2]=chr21Ethan_repC.RNA
FILENAMES[3]=chr21Eric_repA.RNA
FILENAMES[4]=chr21Eric_repB.RNA
FILENAMES[5]=chr21Eric_repC.RNA

for index in $(seq 0 5)
do
rootname=${FILENAMES[$index]}

maindir=/Users/maallen3/sread2021/day5/
indir=${maindir}fastq/
fastqcoutdir=${maindir}fastqc/
trimoutdir=${maindir}trimmomatic/
hisat2outdir=${maindir}hisat2/
genomefile=/scratch/Shares/public/sread2021/data_files/genome/hg38/HISAT2/genome
filename1=${trimoutdir}${rootname}output_forward_paired.fq.gz
filename2=${trimoutdir}${rootname}output_reverse_paired.fq.gz
bamfiledir=${maindir}bams/
sortedbamfiledir=${maindir}sortedbams/

mkdir -p $bamfiledir
mkdir -p $sortedbamfiledir

#/Users/maallen3/sread2021/day5/fastq/chr21Eric_repA.RNA.end1.fastq
#/Users/maallen3/sread2021/day5/trimmomatic/chr21Eric_repA.RNAoutput_forward_paired.fq.gz


fastqc -o $fastqcoutdir ${indir}${rootname}.end1.fastq
fastqc -o $fastqcoutdir ${indir}${rootname}.end2.fastq

java -jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar PE ../fastq/${rootname}.end1.fastq ../fastq/${rootname}.end2.fastq ${trimoutdir}${rootname}output_forward_paired.fq.gz ${trimoutdir}${rootname}output_forward_unpaired.fq.gz ${trimoutdir}${rootname}output_reverse_paired.fq.gz ${trimoutdir}${rootname}output_reverse_unpaired.fq.gz ILLUMINACLIP:/opt/trimmomatic/0.36/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

fastqc -o $fastqcoutdir ${filename1} 
fastqc -o $fastqcoutdir ${filename2}


hisat2 -p 4 -x $genomefile -1 ${filename1} -2 ${filename2} >${hisat2outdir}${rootname}.sam 2>${hisat2outdir}${rootname}.sam.stderr

samtools view -@ 8 -b -o ${bamfiledir}${rootname}.bam ${hisat2outdir}${rootname}.sam
samtools sort -@ 8 -o ${sortedbamfiledir}${rootname}.sorted.bam ${bamfiledir}${rootname}.bam
samtools index ${sortedbamfiledir}${rootname}.sorted.bam ${sortedbamfiledir}${rootname}.sorted.bam.bai


done        
