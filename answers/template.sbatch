#!/bin/bash
#SBATCH --job-name=<JOB-NAME>                     # Job name
#SBATCH --mail-type=ALL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<EMAIL>                       # Where to send mail
#SBATCH --nodes=<n>                               # Number of nodes requested
#SBATCH --ntasks=<n>                              # Number of CPUs (processor cores/tasks)
#SBATCH --mem=<n>gb                               # Memory limit
#SBATCH --time=<00:00:00>                         # Time limit hrs:min:sec
#SBATCH --partition=compute                       # Partition/queue requested on server
#SBATCH --output=/scratch/Users/<USERNAME>/eofiles/<JOB-NAME>.%j.out       # Standard output
#SBATCH --error=/scratch/Users/<USERNAME>/eofiles/<JOB-NAME>.%j.err        # Standard error log

### Displays the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

### Assigns path variables 
INPUT_DIRECTORY=<PATH_TO_INPUT_FILE>
OUTPUT_DIRECTORY=<PATH_TO_OUTPUT_FILE>

### Loads modules
<MODULES_TO_LOAD>

### <SOFTWARE SPECIFICS>

echo Job finished at `date +"%T %a %d %b %Y"`
