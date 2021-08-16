#!/bin/bash
#SBATCH --job-name=<JOB_NAME>                           # Job name
#SBATCH --mail-type=FAIL                                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL>              		# Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:15:00                                 # Time limit hrs:min:sec
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH --output=/scratch/Users/<USERNAME>/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<USERNAME>/eofiles/%x_%j.err


################ SET HEADER ##################################################
# Specify --ntasks, --time, --mem, --output, --error


################ SET REQUIRED VARIABLES ##################################################

INDIR=
OUTDIR=

FILENAME=

##########################################################################################
############### LOAD REQUIRED MODULES ####################################################


##########################################################################################
################## PRINT JOB INFO ########################################################

printf "\nfastq Directory: $INDIR"
printf "\nOutput Directory: $OUTDIR"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

##########################################################################################
################## RUN JOB ###############################################################
mkdir -p ${OUTDIR}

cellranger count --id=<id_name_you_prefer>\
                 --transcriptome=/path_to_ref_genome/ \
                 --fastqs=/path_to_fastqs/\
                 --localcores=8 \
                 --localmem=16 \
                 --sample=<sample_name> \

# transcriptome and fastqs only need to point to directory
# If [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz is your fastq file, --sample=Sample Name.
# cell ranger will find all the files with the same sample name.
# No space between '=' and the part you have to fillout, or it will fail.
                                                                                
echo Job finished at `date +"%T %a %d %b %Y"`
