#!/bin/bash
#SBATCH --job-name=HMMRATAC # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=23:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/maallen3/e_and_o/HMMRATAC.%j.out # Standard output
#SBATCH --error=/scratch/Users/maallen3/e_and_o/HMMRATAC.%j.err # Standard error log


#YOu need to use samtools/1.8, bedtools/2.25.0, python/2.7.14/MACS/2.1.1, python/2.7.14/pandas/0.18.1
#load them here

#you will need to run HMMRATAC so you will need the path for that program
HMMRATACpath=/scratch/Shares/public/sread2021/algorithms/HMMRATAC/

#maybe you should set indir, outdir, and rootname as a variable?
#a bam filename variable might be nice
#and did you make the outdirectory yet?

#Your going to need a .genome.info file. The https://github.com/LiuLabUB/HMMRATAC/ has instructions on making that file

#run it using java here. Again https://github.com/LiuLabUB/HMMRATAC/ can help yo. n
java -jar ${HMMRATACpath}HMMRATAC_V1.2.10_exe.jar -b ${indir}${rootname}.sorted.bam -i ${indir}${rootname}.sorted.bam.bai -g ${outdir}${rootname}.genome.info -o ${outdir}${rootname}

#mabe you should filer the peaks with awk like they suggested https://github.com/LiuLabUB/HMMRATAC/

#crap, are the peaks in order?

