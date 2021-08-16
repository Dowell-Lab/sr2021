#!/bin/bash
#SBATCH --job-name=<JOB_NAME>                           # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL>                        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:30:00                                 # Time limit hrs:min:sec
#SBATCH --partition=compute                             # Partition/queue requested on server    
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH --output=/scratch/Users/<USERNAME>/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<USERNAME>/eofiles/%x_%j.err


################### LOAD MODULES ###################################

module load python/2.7.14
module load python/2.7.14/MACS/2.1.1
module load bedtools/2.25.0

################## SET VARIABLES ###################################
###Joes Variables
ControlBAM='/scratch/Workshop/SR2019/ChIP-seq/mapped/bams/SRR5855054_chr1.sorted.bam'
BAMS_DIRECTORY='/scratch/Workshop/SR2019/ChIP-seq/mapped/bams/'
###COND1='SRR5855055_chr1'
COND1='SRR5855056_chr1'
OUT='/scratch/Users/joefranchesco/ChIP/macs2'
BLACKLIST='/scratch/Workshop/hg38/wgEncodeDacMapabilityConsensusExcludable_grch38.bed'
ROOTNAME=$COND1


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

#######################################################################

echo macs2
date
date

macs2 callpeak \
 -c ${ControlBAM} \
 -t $BAMS_DIRECTORY/${COND1}.sorted.bam \
 --outdir ${OUT} \
 -n ${ROOTNAME} \
 -g hs \
 --SPMR \
 --broad \
 -q 0.00001


echo removing blacklist regions
date
date

bedtools intersect \
 -a $OUT/${ROOTNAME}_peaks.broadPeak \
 -b ${BLACKLIST} \
 -v \
 > $OUT/${ROOTNAME}_peaks_clean.broadPeak

echo macs2 done
date
date


