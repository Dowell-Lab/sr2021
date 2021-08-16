#!/bin/bash
#SBATCH --job-name=<JOB_NAME>                           # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL>              		# Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=4                                      # Number of CPU (tasks)
#SBATCH --time=12:00:00                                 # Time limit hrs:min:sec
#SBATCH --partition=short                               # Partition/queue requested on server    
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH --output=/scratch/Users/<USERNAME>/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<USERNAME>/eofiles/%x_%j.err


##########################################################################################
############### LOAD REQUIRED MODULES ####################################################
#we will want to specify software requirements

module purge
module load samtools/1.8    # for working with bam files
module load bedtools/2.25.0 # generating coverage files (i.e. bedgraphs)
module load gcc/7.1.0       # required for running Tfit
module load openmpi/1.6.4   # required for running Tfit


################ SET REQUIRED VARIABLES ##################################################
## the input files for this script are bam files
## additionally, we will also want to specify paths to executables

#Tfit paths
TFIT=/path/to/Tfit/executables
TFIT_CONFIG=/path/to/Tfit/config_file

#Input and output directories
BAM=/path/to/bams
BEDGRAPH=/path/to/bedgraphs
TFIT_OUT=/where/to/save/tfit/output

#genome files
CHROM_SIZES=/path/to/genome.chrom.sizes

FILENAME=dmso_rep1.chr21

###########################################################
################## PRINT JOB INFO #########################

printf "\nfastq Directory: $INDIR"
printf "\nOutput Directory: $OUTDIR"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

########################################################
#########            RUN JOB               #############
########################################################


###########################################                                                                                                      
##1: Generate bam files and filter chr21 ##
###########################################

samtools view -@ 4 -h -q 1 ${BAM}/${FILENAME}.bam \
    | awk 'substr($1,0,1)=="@" || $0 !~ /ZS:/' \
    | samtools view -h -b \
    > ${BAM}/${FILENAME}.filtered.bam

samtools index ${BAM}/${FILENAME}.filtere.bam \
    ${BAM}/${FILENAME}.filtered.bam.bai


#########################################
##2 : Generate bedgraphs for TFit      ##
#########################################

genomeCoverageBed \
    -bg \
    -strand + \
    -g ${CHROM_SIZES} \
    -ibam ${BAM}/${FILENAME}.filtered.bam \
    > ${BEDGRAPH}/${FILENAME}.pos.bedGraph

genomeCoverageBed \
    -bg \
    -strand - \
    -g ${CHROM_SIZES} \
    -ibam ${BAM}/${FILENAME}.filtered.bam \
    | awk 'BEGIN{FS=OFS="\t"} {$4=-$4}1' \
    > ${BEDGRAPH}/${FILENAME}.neg.bedGraph

cat ${BEDGRAPH}/${sample}.pos.bedGraph \
    ${BEDGRAPH}/${FILENAME}.neg.bedGraph \
    > ${BEDGRAPH}/${FILENAME}.unsorted.bedGraph

sortBed \
    -i ${BEDGRAPH}/${FILENAME}.unsorted.bedGraph \
    > ${BEDGRAPH}/${FILENAME}.bedGraph


#########################################
##3: Run Tfit                         ###
######################################### 

mkdir -p ${TFIT_OUT}
mpirun -np 1 -host ${SLURM_JOB_NODELIST} \
    ${TFIT} bidir -config ${TFIT_CONFIG} \
    -ij ${BEDGRAPH}/${FILENAME}.bedGraph \
    -N ${FILENAME} -o ${TFIT_OUT}

########################################################
#########             DONE!                #############
########################################################

                                                                                
echo Job finished at `date +"%T %a %d %b %Y"`
