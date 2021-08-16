#!/bin/bash
#SBATCH --job-name=fastqc                     # Job name
#SBATCH --mail-type=ALL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<your_email>                    # Where to send mail
#SBATCH --nodes=1                               # Number of nodes requested
#SBATCH --ntasks=1                              # Number of CPUs (processor cores/tasks)
#SBATCH --mem=10gb                               # Memory limit
#SBATCH --time=00:20:00                        # Time limit hrs:min:sec
#SBATCH --partition=short                       # Partition/queue requested on server
#SBATCH --output=</path/to/your/error_directory>%x.%j.out       # Standard output
#SBATCH --error=</path/to/your/error_directory>%x.%j.err        # Standard error log

### Displays the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

### Assigns path variables 
INDIR=</path/to/fastq_dir>
OUTDIR=</path/to/your/FastQC_dir>

### Loads modules
module load fastqc

###Check FASTQ Files
mkdir -p ${OUTDIR}/fastqc
fastqc ${INDIR}/<file.fastq> -o ${OUTDIR}/fastqc

echo Job finished at `date +"%T %a %d %b %Y"`
