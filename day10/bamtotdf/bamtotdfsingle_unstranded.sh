#!/bin/bash 
#SBATCH --job-name=bamtotdfsingle_unstranded # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=username@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=1gb # Memory limit
#SBATCH --time=01:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/<username>/e_and_o/bamtotdf.%j.out # Standard output
#SBATCH --error=/scratch/Users/<username>/e_and_o/bamtotdf.%j.err # Standard error log

module load python_2.7.3
module load samtools_0.1.19
module load bedtools2_2.22.0
module load igvtools_2.1.24



bamroot=#rootfilename
id=#indir
od=#outdir
bgd=#bedgraphdir
bgft=#bedgraphfortdfdir
bgg=#bedgraphgenomefile
bd=#scriptdirname
tdfd=#tdfdir
igvgenome=#igvtoolsgenomefile



samtools sort -m500000000000 ${id}${bamroot}.bam ${od}${bamroot}.sorted
echo finished bam sort
samtools flagstat ${id}${bamroot}.bam > ${od}${bamroot}.bam.flagstat 2>${od}${bamroot}.bam.flagstat.err
echo finished counting infile
samtools index ${od}${bamroot}.sorted.bam
echo finished indexing sorted bam file
samtools/0.1.19/samtools flagstat ${od}${bamroot}.sorted.bam > ${od}${bamroot}.sorted.bam.flagstat 2>${od}${bamroot}.sorted.bam.flagstat.err
echo finished counting sorted bam file
diff ${outdir}${rootfilename}.sorted.bam.flagstat ${outdir}${rootfilename}.bam.flagstat
echo finished checking that the sorted bam file has the same number of reads as the unsorted
bedtools genomeCoverageBed -bg -strand + -ibam ${od}${bamroot}.sorted.bam -g $bgg > ${bgd}${bamroot}.pos.BedGraph
echo finished building histogram of positive strand
bedtools genomeCoverageBed -bg -strand - -ibam ${od}${bamroot}.sorted.bam -g $bgg | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }'> ${bgd}${bamroot}.neg.BedGraph
echo finished building histogram of negative strand
cat ${bgd}${bamroot}.pos.BedGraph ${bgd}${bamroot}.neg.BedGraph > ${bgft}${bamroot}.BedGraph.temp
echo finished concatining Bedgraph files
bedtools sortBed -i ${bgft}${bamroot}.BedGraph.temp >${bgft}${bamroot}.BedGraph
echo finished sorting the Begraph file with both strands
python ${bd}readcount_corrected_geneomeBedgraphs.py ${od}${bamroot}.sorted.bam.flagstat ${bgft}${bamroot}.BedGraph
echo finished correcting for millions mapped
igvtools toTDF ${bgft}${bamroot}.BedGraph.mp.BedGraph ${tdfd}${bamroot}.tdf $igvgenome
echo finished making tdf



