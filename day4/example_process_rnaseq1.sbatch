#!/bin/bash
#SBATCH --job-name=example_rnaseq_process       # Job name
#SBATCH --mail-type=ALL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL_HERE>           # Where to send mail
#SBATCH --nodes=1                               # Number of nodes requested
#SBATCH --ntasks=8                              # Number of CPUs (processor cores/tasks)
#SBATCH --mem=2gb                               # Memory limit
#SBATCH --time=01:30:00                        # Time limit hrs:min:sec
#SBATCH --partition=short                       # Partition/queue requested on server
#SBATCH --output=/scratch/Shares/dowell/sread/scripts/day4/e_and_o/%x.%j.out       # Standard output:  the %x will be replaced with job_name and the %j by the job id
#SBATCH --error=/scratch/Shares/dowell/sread/scripts/day4/e_and_o/%x.%j.err        # Standard error log

### Displays the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

##############################
### Assigns path variables ###
##############################

## the fastq files will be used as input to fastqc and trimmomatic
## trimmed reads will then be passed on to the mapping step
DATADIR=/scratch/Shares/dowell/sread/data_files/day4
FASTQ=${DATADIR}/fastq
FASTQC=${DATADIR}/fastqc
HISAT2=${DATADIR}/hisat2
TRIM=${DATADIR}/trimmomatic

FILENAME=chr21Eric_repA
######################
### Loads modules   ##
######################
module load fastqc 
module load trimmomatic
module load hisat2/2.1.0
module load samtools/1.8

#######################
##Run the pipeline   ##
#######################
##1: Run fastqc on the samples (here run on example file ${FILENAME}.RNA.end1.fastq)
fastqc ${FASTQ}/${FILENAME}.RNA.end1.fastq -o ${FASTQC}
fastqc ${FASTQ}/${FILENAME}.RNA.end2.fastq -o ${FASTQC}

##2: Trim FASTQ Files                                                                                                                                                                                                                                                
java -jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar PE \
 -threads 8 \
 -phred33 \
 -trimlog ${TRIM}/trimlog \
 ${FASTQ}/${FILENAME}.RNA.end1.fastq ${FASTQ}/${FILENAME}.RNA.end2.fastq \
 ${TRIM}/${FILENAME}.RNA.end1.trimmed.fastq ${TRIM}/${FILENAME}.RNA.end1.unpaired.fastq \
 ${TRIM}/${FILENAME}.RNA.end2.trimmed.fastq ${TRIM}/${FILENAME}.RNA.end2.unpaired.fastq \
 ILLUMINACLIP:/opt/trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:30:10 \
 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

##3: Check Post-Trimming QC stats
fastqc ${TRIM}/*.trimmed.fastq -o ${FASTQC}

##4: Map trimmed reads to reference genome                                                                                                                                                     
hisat2 --very-fast -x /scratch/Shares/public/genomes/hisatfiles/hg38/HISAT2/genome \
-1 ${TRIM}/${FILENAME}.RNA.end1.trimmed.fastq \
-2 ${TRIM}/${FILENAME}.RNA.end2.trimmed.fastq \
> ${HISAT2}/${FILENAME}.RNA.sam \
2> ${HISAT2}/${FILENAME}.hisat2_maptstats.txt

##5: Convert SAM file to BAM file and index the BAM files
##The BAM file is a compressed binary version of a SAM file, so takes up less space
##The sorted bam also takes up less space
###5.1 - convert SAM to BAM 
samtools view -@ 8 -bS -o ${HISAT2}/${FILENAME}.RNA.bam ${HISAT2}/${FILENAME}.RNA.sam

###5.2 - sort bam file
samtools sort -@ 8 ${HISAT2}/${FILENAME}.RNA.bam > ${HISAT2}/${FILENAME}.RNA.sorted.bam

###5.3 - index sorted bam file
samtools index ${HISAT2}/${FILENAME}.RNA.sorted.bam  ${HISAT2}/${FILENAME}.RNA.sorted.bam.bai


#TO DO: In the actual script we can delete the intermediate SAM, and unsorted BAM file
 
echo Job finished at `date +"%T %a %d %b %Y"`
