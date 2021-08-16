#!/bin/bash
#SBATCH --job-name=featurecounts                                 # Job name
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu                 # Where to send mail
#SBATCH --nodes=1                                        # Number of cores job will run on
#SBATCH --ntasks=4                                       # Number of CPU (processers, tasks)
#SBATCH --time=2:00:00                                   # Time limit hrs:min:sec
#SBATCH --partition compute                              # Job queue
#SBATCH --mem=4gb                                        # Memory limit
#SBATCH --output=/scratch/Users/maallen3/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/maallen3/eofiles/%x_%j.err



module load R/3.6.1
export R_LIBS="/data/R-lib" #this is here becuase we don't all want to install feature counts on R in our own home directorys. BIT installed it here. 

################### SET VARIABLES ######################################

FEATURECOUNTS=featureCounts.R

########################################################################
################## PRINT JOB INFO ######################################

printf "Sample ID: $ROOTNAME"
printf "\nDirectory: $PROJECT"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

#######################################################################

Rscript $FEATURECOUNTS
