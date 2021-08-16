#!/bin/bash 
#SBATCH --job-name=bamtotdfpaired # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=8     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=01:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/maallen3/e_and_o/bamtotdfpaired.%j.out # Standard output
#SBATCH --error=/scratch/Users/maallen3/e_and_o/bamtotdfpaired.%j.err # Standard error log


module load samtools/1.8
module load bedtools/2.25.0
module load python/2.7.14/
module load igvtools/2.3.75

bamroot=SRR11856163
id=/scratch/Shares/public/sread2021/data_files/day10/filetypes/
od=/scratch/Users/maallen3/day10/bamtotdf/single/
genomefasta=/scratch/Shares/public/genomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
bgft=${od}bedgraphfortdf/
tdfd=${od}tdfdir/
statdir=${od}stats/
mkdir -p $bgft $tdfd $statdir


#get the chromosome sizes files you will need to do this
samtools faidx $genomefasta | cut -f1,2 > ${statdir}chrom.sizes
igvgenome=${statdir}chrom.sizes
echo got the igvgenome and bedtools genomefile


mainbam=${id}${bamroot}.sorted.bam
samtools flagstat -@ 8 $mainbam > ${statdir}${bamroot}.bam.flagstat 2>${statdir}${bamroot}.bam.flagstat.err
echo finished counting infile
flagstatfile=${statdir}${bamroot}.sorted.bam.flagstat
echo finished counting sorted bam file


bedtools genomeCoverageBed -bg -strand + -ibam $mainbam -g $igvgenome > ${bgft}${bamroot}.pos.BedGraph
echo finished building histogram of positive strand
bedtools genomeCoverageBed -bg -strand - -ibam $mainbam -g $igvgenome | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }'> ${bgft}${bamroot}.neg.BedGraph
echo finished building histogram of negative strand
cat ${bgft}${bamroot}.pos.BedGraph ${bgft}${bamroot}.neg.BedGraph > ${bgft}${bamroot}.BedGraph.temp
echo finished concatining Bedgraph files


bedtools sortBed -i ${bgft}${bamroot}.BedGraph.temp >${bgft}${bamroot}.BedGraph
echo finished sorting the Begraph file with both strands
python ${bd}readcount_corrected_geneomeBedgraphs.py $flagstatfile ${bgft}${bamroot}.BedGraph
echo finished correcting for millions mapped
igvtools toTDF ${bgft}${bamroot}.BedGraph.mp.BedGraph ${tdfd}${bamroot}.tdf $igvgenome
echo finished making tdf


