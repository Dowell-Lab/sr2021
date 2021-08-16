#!/bin/bash
#SBATCH --job-name=hisat2_mapping                 # Job name
#SBATCH --mail-type=ALL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=qing.yang@colorado.edu        # Where to send mail
#SBATCH --nodes=1                                 # Number of nodes requested
#SBATCH --ntasks=8                                # Number of CPUs (processor cores/tasks)
#SBATCH --mem=2gb                                 # Memory limit
#SBATCH --time=01:30:00                         # Time limit hrs:min:sec
#SBATCH --partition=short                       # Partition/queue requested on server
#SBATCH --output=/scratch/Users/qiya9811/eofiles/test.%j.out       # Standard output
#SBATCH --error=/scratch/Users/qiya9811/eofiles/test.%j.err        # Standard error log

### Displays the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

### Assigns path variables

DATADIR=/scratch/Shares/dowell/sread/cookingShow/day4/
HISAT2=${DATADIR}/hisat2
TRIM=${DATADIR}/trimmomatic

### Loads modules
module load hisat2/2.1.0
module load samtools/1.8

### <SOFTWARE SPECIFICS>

## Map trimmed reads to reference genome                                                                                                                                                     
hisat2 --very-fast -x /scratch/Shares/public/genomes/hisatfiles/hg38/HISAT2/genome \
-1 ${TRIM}/chr21Eric_repA.RNA.end1.trimmed.fastq \
-2 ${TRIM}/chr21Eric_repA.RNA.end2.trimmed.fastq \
> ${HISAT2}/chr21Eric_repA.RNA.sam \
2> ${HISAT2}/chr21Eric_repA.hisat2_maptstats.txt

## Convert mapped reads to sorted bam file
### convert SAM to BAM
samtools view -@ 8 -bS -o ${HISAT2}/chr21Eric_repA.RNA.bam \
${HISAT2}/chr21Eric_repA.RNA.sam

### sort bam file
samtools sort -@ 8 ${HISAT2}/chr21Eric_repA.RNA.bam \
> ${HISAT2}/chr21Eric_repA.RNA.sorted.bam

### index sorted bam file
samtools index ${HISAT2}/chr21Eric_repA.RNA.sorted.bam \
${HISAT2}/chr21Eric_repA.RNA.sorted.bam.bai


echo Job finished at `date +"%T %a %d %b %Y"`
