#!/bin/bash 
#SBATCH -p compute
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=username@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node 
#SBATCH --ntasks=4     # Number of CPU (processer cores i.e. tasks)
#SBATCH --mem=5gb # Memory limit
#SBATCH --time=00:15:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/username/e_and_o/dmso_rep1.chr1_Tfit.%j.out
#SBATCH --error=/scratch/Users/username/e_and_o/dmso_rep1.chr1_Tfit.%j.err
#SBATCH --job-name=Tfit_run # Job name

#load modules
module load mpi/openmpi-x86_64
module load gcc/7.1.0
module load bedtools/2.25.0

#Initiate paths and variables
outdir=/scratch/Users/username/day10
rootname=dmso_rep1.chr21
bdgraphdir=/scratch/Shares/public/sread2021/cookingShow/day8/bedgraph_groseq

###########################################################                        
################## print job info #########################                                                                                                                                                     
printf "\nfastq Directory: $INDIR"
printf "\nOutput Directory: $OUTDIR"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"
date

export OMP_NUM_THREADS=4
echo $OMP_NUM_THREADS
echo $SLURM_JOB_NODELIST 

#make output directories
mkdir -p $outdir
mkdir -p ${outdir}/tfit_out

tfit_outdir=${outdir}/tfit_out
bg_file=${bdgraphdir}/${rootname}.bedGraph

#Initiate Tfit executables
Tfitdir=/scratch/Shares/public/sread2021/algorithms/Tfit_2017
src=${Tfitdir}/src/Tfit
config_file=${Tfitdir}/config_files/config_file.txt

#================================================================
#calling tfit command
cmd="mpirun -np 1 -host ${SLURM_JOB_NODELIST}"

$cmd $src bidir -ij $bg_file -o $tfit_outdir -N ${rootname} -MLE 1
#================================================================
