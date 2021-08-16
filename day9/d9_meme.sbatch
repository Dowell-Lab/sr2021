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

module load bedtools/2.25.0
module load meme/5.1.1

################## SET VARIABLES ###################################

#INDIR is where your macs output peak files are stored. OUTDIR is where I want my output from running meme and tomtom
INDIR='/path/to/macs/peak/files'
OUTDIR='/path/to/output'
mkdir -p ${OUTDIR}

HG38_FASTA='/scratch/Shares/dowell/genomes/hg38/hg38.fa'

FILENAME='BACH1'

################## JOB INFO ########################################

printf "\nMACS Peak Output Directory: ${INDIR}"
printf "\nOutput Directory: ${OUTDIR}"
printf "\nOutput File Name: $FILENAME"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

################### RUN JOB ########################################

#### Get fasta of peak files
echo convert peaks call to fasta format
date 
date

bedtools getfasta \
 -fi ${HG38_FASTA} \
 -bed ${OUTDIR}/${FILENAME}_peaks_clean.narrowPeak \
 -fo ${OUTDIR}/${FILENAME}_peaks_clean.fasta

echo fasta file done 
date
date

#### Run MEME command line
# URL https://meme-suite.org/meme/tools/meme
#/opt/meme/5.1.1/bin/meme \
# ${INDIR}/${ROOTNAME}.fasta \
# -dna \
# -oc ${OUTDIR}/meme/${ROOTNAME} \
# -nostatus \
# -time 18000 \
# -mod zoops \
# -nmotifs 10 \
# -minw 6 \
# -maxw 50 \
# -objfun classic \
# -revcomp \
# -markov_order 0


### Run TOMTOM command line
# URL https://meme-suite.org/meme/tools/tomtom
#/opt/meme/5.1.1/bin/tomtom \
# ${OUTDIR}/meme/${ROOTNAME}/meme.html \
# -no-ssc \
# -oc ${OUTDIR}/tomtom/${ROOTNAME} \
# -verbosity 1 \
# -min-overlap 5 \
# -dist pearson \
# -evalue \
# -thresh 10.0
