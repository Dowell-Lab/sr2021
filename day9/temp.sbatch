#!/bin/bash
#SBATCH --job-name=macs2                               # Job name
#SBATCH --mail-type=ALL                                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<USER@EMAIL.COM>                   # Where to send mail
#SBATCH --nodes=1                                      # Numbers of nodes
#SBATCH --ntasks=1                                     # Number of CPU (tasks)
#SBATCH --time=00:30:00                                # Time limit hrs:min:sec
#SBATCH --partition=compute                            # Partition/queue requested on server    
#SBATCH --mem=2gb                                      # Memory limit
#SBATCH --output=/path/to/eofiles/%x_%j.out
#SBATCH --error=/path/to/eofiles/%x_%j.err


################### LOAD MODULES ###################################

module load python/2.7.14
module load python/2.7.14/MACS/2.1.1
module load bedtools/2.25.0

################## SET VARIABLES ###################################

#INDIR is where my bams file are stored. OUTDIR is where I want my output from running MACS to go
INDIR='/path/to/bam/file'
OUTDIR='/path/to/output/diretory'

BLACKLIST='/scratch/Shares/dowell/sread/data_files/day9/blacklist_hg38/problematic_regions_hg38.bed'
HG38_FASTA='/scratch/Shares/public/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'

FILENAME='BACH1'


################## JOB INFO ########################################

printf "\nMACS Directory: $INDIR"
printf "\nOutput Directory: $OUTDIR"
printf "\nOutput Rootname: $FILENAME"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

################### RUN JOB ########################################
#Run MEME command line # URL https://meme-suite.org/meme/tools/meme
/opt/meme/5.1.1/bin/meme \
 ${INDIR}/${ROOTNAME}.fasta \
 -dna \
 -oc ${OUTDIR}/meme/${ROOTNAME} \
 -nostatus \
 -time 18000 \
 -mod zoops \
 -nmotifs 10 \
 -minw 6 \
 -maxw 50 \
 -objfun classic \
 -revcomp \
 -markov_order 0


### Run TOMTOM command line
# URL https://meme-suite.org/meme/tools/tomtom
/opt/meme/5.1.1/bin/tomtom \
 ${OUTDIR}/meme/${ROOTNAME}/meme.html \
 -no-ssc \

