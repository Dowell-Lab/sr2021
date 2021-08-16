#!/bin/bash 
#SBATCH --job-name=slurm_test # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=username@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=1gb # Memory limit
#SBATCH --time=01:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/<username>/e_and_o/slurm_test.%j.out # Standard output
#SBATCH --error=/scratch/Users/<username>/e_and_o/slurm_test.%j.err # Standard error log


#The first line in this script tells the computer you want bash (bourne again shell)to run this job.
#the top is the SBATCH information
#the bottom is the command

pwd; hostname; date

echo "You've requested $SLURM_CPUS_ON_NODE core(s)."
sleep 600
date


