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
od=/scratch/Users/maallen3/day10/bamtotdf/
genomefasta=/scratch/Shares/public/genomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
bgft=${od}bedgraphfortdf/
tdfd=${od}tdfdir/
statdir=${od}stats/
mkdir -p $bgft $tdfd $statdir



#Assumption 1, we are starting with sorted indexed bam files
#Assumption 2, Read1 is on the incorrect strand, Read2 is on the correct strand. This is true for most of the kits we use.


mainbam=${id}${bamroot}.sorted.bam
samtools flagstat -@ 8 $mainbam > ${statdir}${bamroot}.bam.flagstat 2>${statdir}${bamroot}.bam.flagstat.err
echo finished counting infile
flagstatfile=${indir}${bamroot}.sorted.bam.flagstat


#get the chromosome sizes files you will need to do this
samtools faidx $genomefasta | cut -f1,2 > ${statdir}chrom.sizes
igvgenome=${statdir}chrom.sizes
echo Got chromosome size
#count millions mapped. Flagstat reports to many reads as mapping when you use it on paired end data. 
samtools view -cF 0x100 $mainbam >${statdir}${bamroot}.bam.millionsmapped
echo finished counting sorted bam with flags

#http://davetang.org/wiki/tiki-index.php?page=SAMTools#Extracting_only_the_first_read_from_paired_end_BAM_files
samtools view -h -b -f 0x0040 ${mainbam} > ${bgft}${bamroot}.pairfirst.bam
# 0x0040 is hexadecimal for 64 (i.e. 16 * 4), which is binary for 1000000, corresponding to the read in the first read pair.
echo finished pulling out the first read with flags
#need to know the flag for the second strand
# Jess gave me this https://broadinstitute.github.io/picard/explain-flags.html
#128 means second in pair
#128 in hexadecimal is 0x0080
samtools view -h -b -f 0x0080 ${mainbam} > ${bgft}${bamroot}.pairsecond.bam
echo finished pulling out the second read with flags
#That should get me the two parts of pair separate
#Then I need to run genomecoverage and swap the strand info on the second pair


#Then I need to run genomecoverage on each of them
genomeCoverageBed -bg -split -strand - -ibam ${bgft}${bamroot}.pairfirst.bam -g $igvgenome >${bgft}${bamroot}.pairfirst.pos.bed
echo putting the read1 negative strand reads on the postive strand in the bedgraph
genomeCoverageBed -bg  -split -strand + -ibam ${bgft}${bamroot}.pairfirst.bam -g $igvgenome | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }' >${bgft}${bamroot}.pairfirst.neg.bed
echo putting the read1 postive strand reads on the negative strand in the bedgraph
genomeCoverageBed -bg -split -strand + -ibam  ${bgft}${bamroot}.pairsecond.bam -g $igvgenome > ${bgft}${bamroot}.pairsecond.pos.bed
echo putting the read2 postive strand reads on the postive strand in the bedgraph
genomeCoverageBed -bg -split -strand - -ibam ${bgft}${bamroot}.pairsecond.bam -g $igvgenome | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }'> ${bgft}${bamroot}.pairsecond.neg.bed
echo putting the read2 negative strand reads on the negative strand in the bedgraph
#first I need to sort the Bedgraphs
sortBed -i ${bgft}${bamroot}.pairfirst.pos.bed > ${bgft}${bamroot}.pairfirst.pos.BedGraph.sort

sortBed -i ${bgft}${bamroot}.pairfirst.neg.bed > ${bgft}${bamroot}.pairfirst.neg.BedGraph.sort


sortBed -i ${bgft}${bamroot}.pairsecond.pos.bed > ${bgft}${bamroot}.pairsecond.pos.BedGraph.sort

sortBed -i ${bgft}${bamroot}.pairsecond.neg.bed > ${bgft}${bamroot}.pairsecond.neg.BedGraph.sort
echo finished sorting the bedgraphs

#Then I need to add the two Bedgraphs 

#this should put the values in columns 4 and 5


unionBedGraphs -i ${bgft}${bamroot}.pairfirst.pos.BedGraph.sort ${bgft}${bamroot}.pairsecond.pos.BedGraph.sort >${bgft}${bamroot}.pos.Bedgraphcol

unionBedGraphs -i ${bgft}${bamroot}.pairfirst.neg.BedGraph.sort ${bgft}${bamroot}.pairsecond.neg.BedGraph.sort >${bgft}${bamroot}.neg.Bedgraphcol

#then I need to sum cols 4 and 5

awk -F '\t' '{OFS="\t"; print $1,$2,$3,$4+$5;}' ${bgft}${bamroot}.pos.Bedgraphcol >${bgft}${bamroot}.pos.Bedgraph

awk -F '\t' '{OFS="\t"; print $1,$2,$3,$4+$5;}' ${bgft}${bamroot}.neg.Bedgraphcol >${bgft}${bamroot}.neg.Bedgraph

#then I need to cat the two Bedgraphs
echo finished adding the postive and negative reads from both strands back together
cat ${bgft}${bamroot}.pos.Bedgraph ${bgft}${bamroot}.neg.Bedgraph >${bgft}${bamroot}.bed


#then I need to sort the final Bedgraph so it can be divided by millions mapped and converted into tdf
sortBed -i ${bgft}${bamroot}.bed >${bgft}${bamroot}.BedGraph
echo finished making the strand corrected sorted bedgraph
#now correct for millionsmapped
python readcount_corrected_geneomeBedgraphs.py ${statdir}${bamroot}.bam.flagstat ${bgft}${bamroot}.BedGraph

echo finished correcting for millions mapped
igvtools toTDF ${bgft}${bamroot}.mp.BedGraph ${tdfd}${bamroot}.tdf $igvgenome
echo finished making tdf

