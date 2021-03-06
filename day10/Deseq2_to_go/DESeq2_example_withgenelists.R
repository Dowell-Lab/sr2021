# For more information see: http://www.bioconductor.org/help/workflows/rnaseqGene/
library("tidyverse")
library("DESeq2")
#set working dir
workdir <- '/scratch/Users/username/day10/Deseq2_to_go/'
dir.create(workdir, showWarnings = FALSE)
setwd(workdir)
getwd()
outdir <- paste(workdir, 'deseqresults', '/', sep='') ##naming our outdir
dir.create(outdir, showWarnings = FALSE) ###creating the directory

counts <- read.csv("/scratch/Shares/public/sread2021/cookingShow/day8/RNAseqextras/counts/featureCounts.txt", row.names=1, sep="\t")
head(counts) #your rowname should be the gene ids. Your colnames should match some column of your metadata table (in this case filetable)
filetable <- read.csv('/scratch/Shares/public/sread2021/cookingShow/day8/RNAseqextras/meta.txt', sep="\t")
head(counts)
filetable$chr21 <- factor(filetable$chr21)
filetable$bamfiles <- paste0(filetable$Run, ".sorted.bam") #making a column the files
filelist<- filetable$bamfiles #creating a vector that is the file list

#Your metatdata columns and your counts rows must be in the same order!!!!!!
counts <- counts %>% select(as.vector(filetable$bamfiles))

# Generate DESeqDataSet from count matrix generated by featureCounts
ddsMat <- DESeqDataSetFromMatrix(countData = counts, colData = filetable, design=~chr21)
dds <- ddsMat


### Run DESeq on the DESeqDataSet object
DEdds <- DESeq(dds)

### output the results for a specified alpha value
alpha_val <- 0.05
comparison <- c("chr21", "Disomic", "Trisomic")
res <- results(DEdds, alpha = alpha_val, contrast = comparison)

res_shrink <- lfcShrink(DEdds, contrast = comparison, res = res)

### MA plot
name <- "MA_tri_vs_ctrl_DEA"
limits <- c(-10, 10)
pdf(paste0(outdir, name, ".pdf"))
maplot <- plotMA(res_shrink, main="Disomic vs Trisomic", alpha=alpha_val, ylim=limits)
dev.off()

### disp plot
name <- "disp_tri_vs_ctrl_DEA"
limits <- c(-10, 10)
pdf(paste0(outdir, name, ".pdf"))
maplot <- plotDispEsts(DEdds, main="Disomic vs Trisomic")
dev.off()

#### sort by sig
res_shrink<- res_shrink[ order( res_shrink$padj ), ]

### Take subset of results that are significant
res_shrink_Sig <- subset(res_shrink, padj < alpha_val)


write.csv(res_shrink, file = paste0(outdir,"all_results.csv"))
write.csv(res_shrink_Sig, file = paste0(outdir,"sig_results.csv"))


#for go and enricher and gsea
res_shrink_expressed <- as.data.frame(res_shrink)
res_shrink_expressed <- res_shrink_expressed[!is.na(res_shrink_expressed$padj),]
write.csv(rownames(res_shrink_expressed), file = paste0(outdir,"backgroundgenes.csv"),row.names = FALSE, col.names = FALSE, quote = FALSE)
write.csv(rownames(res_shrink_Sig), file = paste0(outdir,"siggenes.csv"),row.names = FALSE, col.names = FALSE, quote = FALSE)

rnkdf <- tibble(gene = rownames(res_shrink),
				rnk = -log(res$pvalue) * sign(res$log2FoldChange)) %>%
	arrange(desc(rnk)) %>% drop_na()

## Write out the table without any additional information
write.table(rnkdf, file = paste0(outdir,"deseq_res_for_gsea.rnk"),
			append = FALSE, col.names = FALSE, row.names = FALSE,
			quote = FALSE, sep = "\t")

