#!/bin/bash 
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=1gb # Memory limit
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --job-name=featurecounts                                 # Job name
#SBATCH --mail-user=email@colorado.edu           # Where to send mail
#SBATCH --partition short                                # Job queue
#SBATCH --output=/scratch/Users/username/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/username/eofiles/%x_%j.err

filename=SRR15283267
outdir=

module load sra/2.9.2

fastq-dump -O $outdir -split-3 $filename
