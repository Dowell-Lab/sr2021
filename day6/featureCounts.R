####### Running featureCounts
####### Author: Taylor Jones
####### Here we will learn how to download a package, what metadata table is (and why it is important),
####### and run featureCounts, which counts reads over genes.

#this line is only because we won't install libs on the AWS
.libPaths("/data/R-lib")

##############################################################################################################################
####### We will want to start fresh and clear our environment.
# start by clearing your console. To do this hit Ctrl+l or go to Edit-->Clear Console
# clear your environment and plots by hitting the broom icon in both those cells 
# reset our working directory
workdir <- '/scratch/Users/username/day6/'
setwd(workdir)
getwd()

# Like before: most packages you can install with this syntax: install.packages("PACKAGE"). Such as this package: install.packages("tidyverse") --great package
# HOWEVER, some of the bioinformatic software is located in Bioconductor, including RSubread which contains featureCounts.

### To install un-comment the next two lines and run them:
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")   #this will take a moment to run
#BiocManager::install("RSubread")  #This actually installs RSubread
### remember to hash out the code when run is complete

library("Rsubread") #this actually loads the library into the environment
# ??Rsubread #this calls for help of a package. This is useful if you have no idea where to start
### click on the Help pages search results
# ??Rsubread::featureCounts #we can look specifically at featureCounts and the featureCount flags

####### NOTE: We run featureCounts on a cluster because bam files are LARGE
####### We want to run this on the cluster but using command line R is not always great.
####### An alternatve is to generate a script like this then submit this script to the cluster in an sbatch script.
####### Read in bam file list for featureCounts
#bamdir <- '/scratch/Users/YOU/day4/hisat2/'
bamdir <-'/scratch/Shares/public/sread2021/cookingShow/day6_7/'

#chr21Eric_repA.RNA.bam
filelist <- c(paste0(bamdir, "chr21Eric_repA.RNA.sorted.bam"),
           paste0(bamdir, "chr21Ethan_repA.RNA.sorted.bam"))

outdir <- paste(workdir, 'counts', '/', sep='') ##naming our outdir
dir.create(outdir, showWarnings = FALSE) ###creating the directory

# We use a GTF -- this is necessary for featureCounts and more useful for RNA-seq data as it is cleaner/faster for counting over exons
hg38gtf <- "/scratch/Shares/public/genomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/chr21_genes.gtf"

# Read counting using featureCounts
# Sink will save out stdout -- check all of these settings (e.g. isPairedEnd)!
sink(paste0(outdir, "featureCounts_gene_rnaseq.notes.txt")) #we want to save the output as a file
fc <- featureCounts(files=filelist,
                    annot.ext=hg38gtf,
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="exon",
                    GTF.attrType="gene_id",
                    useMetaFeatures=TRUE,
                    allowMultiOverlap=TRUE,
                    largestOverlap=TRUE,
                    countMultiMappingReads=TRUE,
                    isPairedEnd=TRUE,
                    nthreads=4)  #when you move to a bigger machine change to 8
#These parameters above are for specifying how to count over the gtf file.
#Flag info here: ??Rsubread::featureCounts for specifics on each flag
sink()

### Write results -- we'll keep our accession number (GeneID) and length for filtering genes in DESeq2
write.table(x=data.frame(fc$annotation[,c("GeneID")],
                         fc$counts,stringsAsFactors=FALSE),
            paste0(outdir, file="featureCounts_gene_rnaseq.txt"),
            quote=FALSE,sep="\t",
            row.names=FALSE)



