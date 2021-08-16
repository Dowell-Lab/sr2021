#!/bin/bash
#SBATCH --job-name=MACS                                 # Job name
#SBATCH --mail-type=FAIL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jewe1055@colorado.edu                        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:30:00                                 # Time limit hrs:min:sec
#SBATCH --partition=short                             # Partition/queue requested on server    
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH --output=/scratch/Users/jewe1055/sread/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/jewe1055/sread/eofiles/%x_%j.err


################### LOAD MODULES ###################################

module load python/2.7.14
module load python/2.7.14/MACS/2.1.1
module load bedtools/2.25.0

################## SET VARIABLES ###################################

#INDIR is where my bams file are stored. OUTDIR is where I want my output from running MACS to go
INDIR='/scratch/Shares/dowell/sread/data_files/day9/bams'
OUTDIR='/scratch/Users/jewe1055/sread/macs'
mkdir -p ${OUTDIR}/w_ctrl ${OUTDIR}/no_ctrl 

BLACKLIST='/scratch/Shares/dowell/sread/data_files/day9/blacklist_hg38/problematic_regions_hg38.bed'
HG38_FASTA='/scratch/Shares/dowell/genomes/hg38/hg38.fa'

FILENAME='BACH1'

################## JOB INFO ########################################

printf "\nControl ID: $CTRL1, $CTRL2"
printf "\nCondition ID: $COND1, $COND2"
printf "\nBAM Directory: $BAMS"
printf "\nOutput Directory: $OUT"
printf "\nOutput Rootname: $ROOTNAME"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

################### RUN JOB ########################################

echo macs2
date
date

#### Call peaks without controls
macs2 callpeak \
 -t ${INDIR}/${FILENAME}.chr21.sorted.bam \
 --outdir ${OUTDIR}/no_ctrl \
 -g hs \
 -n ${FILENAME} \
 -g hs \
 -B \
 -q 0.01

#### Call peaks with controls
macs2 callpeak \
 -c ${INDIR}/${FILENAME}.input.chr21.sorted.bam \
 -t ${INDIR}/${FILENAME}.chr21.sorted.bam \
 --outdir ${OUTDIR}/w_ctrl \
 -n ${FILENAME} \
 -g hs \
 -B \
 -q 0.01 \

echo macs2 done
date
date

#### Removing ENCODE Blacklist regions
echo removing blacklist regions
date
date

bedtools intersect \
 -a ${OUTDIR}/w_ctrl/${FILENAME}_peaks.narrowPeak \
 -b ${BLACKLIST} \
 -v \
 > ${OUTDIR}/w_ctrl/${FILENAME}_peaks_clean.narrowPeak

echo blacklist regions removed
date
date


#### Get fasta of peak files
echo convert peaks call to fasta
date 
date

bedtools getfasta \
 -fi ${HG38_FASTA} \
 -bed ${OUTDIR}/w_ctrl/${FILENAME}_peaks_clean.narrowPeak \
 -fo ${OUTDIR}/w_ctrl/${FILENAME}_peaks_clean.fasta

echo fasta file 
date
date

