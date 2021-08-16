#!/bin/bash
#SBATCH --job-name=HMMRATAC # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=23:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/maallen3/HMMRATAC/e_and_o/HMMRATAC.%j.out # Standard output
#SBATCH --error=/scratch/Users/maallen3/HMMRATAC/e_and_o/HMMRATAC.%j.err # Standard error log


module load samtools/1.8
module load bedtools/2.25.0
module load python/2.7.14/MACS/2.1.1
module load python/2.7.14/pandas/0.18.1

HMMRATACpath=/scratch/Shares/public/sread2021/algorithms/HMMRATAC/

mkdir -p $outdir

samtools view -H ${indir}${rootname}.sorted.bam| perl -ne 'if(/^@SQ.*?SN:(.+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > ${outdir}${rootname}.genome.info

java -jar ${HMMRATACpath}HMMRATAC_V1.2.10_exe.jar -b ${indir}${rootname}.sorted.bam -i ${indir}${rootname}.sorted.bam.bai -g ${outdir}${rootname}.genome.info -o ${outdir}${rootname}

awk -v OFS="\t" '$13>=10 {print}' ${outdir}${rootname}_peaks.gappedPeak > ${outdir}${rootname}.filteredPeaks.gappedPeak
awk -v OFS="\t" '$5>=10 {print}' ${outdir}${rootname}_summits.bed > ${outdir}${rootname}.filteredSummits.bed

python ${HMMRATACpath}HMMRARACopen.py ${outdir}${rootname}_peaks.gappedPeak ${outdir}${rootname}_peaks.open.bed

bedtools sort -i ${outdir}${rootname}_peaks.open.bed > ${outdir}${rootname}.open.sorted.bed #use these for TFEA
bedtools sort -i ${outdir}${rootname}_summits.bed > ${outdir}${rootname}.summits.sorted.bed

