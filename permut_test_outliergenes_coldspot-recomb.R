#### permutation test using regioneR package to test for an enrichment of outlier genes (gene-based baypass output) in coldspots of recombination identified by Morgan et al 2017 Genetics
### Date 07-05-2020
# Author: C. Smadja

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("regioneR")

library(readr)
library(regioneR)

# install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("biocLite"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

library(readr)
library(regioneR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# load coldspot of recombination info and outlier genes info
coldspots <- toGRanges(A="coldspots.tsv")

########### 0.05 outliers 
outliergenes <- read.table("musCNC_final_outliers_0.05_table_outlierinfo_noheader_chr_sorted.bed", sep = "\t", header = FALSE)
colnames(outliergenes) <- c("chr","start","end","geneID","geneName","significance_criteria", "length","SNP_number","pvalCmax","pvalCmean","status")


### 10000 permutations 06/05/20 ####

# permutation test
chromosomes = outliergenes$chr
starts <- outliergenes$start
ends   = outliergenes$end
minoverlap = 1
maxTrials = 10000

widths = ends - starts

meanWidth = mean(widths)
nbCandidates = length(starts)
df  = data.frame(chrom=chromosomes, starts = starts, ends = ends)

gr  = GRanges(seqnames=chromosomes,ranges=IRanges(start=df$starts, end=df$ends))

numOverlaps(gr, coldspots, count.once=TRUE)
# [1] 186

Recspt <- overlapPermTest(A=gr, B=coldspots, ntimes = 10000, genome="mm10", alternative="auto",force.parallel = FALSE)
pdf(file="permuTest_outliergenes_0.05_coldspots_auto_10000permut.pdf")
plot(Recspt)
dev.off()
lz <- localZScore(A=gr, pt=Recspt, window=10*mean(width(gr)),
                  step=mean(width(gr))/2, B=coldspots)
pdf(file="permuTest_outliergenes_0.05_coldspots_auto_10000permut_Zplot.pdf")
plot(lz)
dev.off()

Recspt <- overlapPermTest(A=gr, B=coldspots, ntimes = 10000, genome="mm10", alternative="greater",force.parallel = FALSE)
# The masked version of 'mm10' is not installed. Using the unmasked version. This means that no automatic masking will be available.
pdf(file="permuTest_outliergenes_0.05_coldspots_greater_10000permut.pdf")
plot(Recspt)
dev.off()
lz <- localZScore(A=gr, pt=Recspt, window=10*mean(width(gr)),
                  step=mean(width(gr))/2, B=coldspots)
pdf(file="permuTest_outliergenes_0.05_coldspots_greater_10000permut_Zplot.pdf")
plot(lz)
dev.off()

Recspt <- overlapPermTest(A=gr, B=coldspots, ntimes = 1000, genome="mm10", alternative="less",force.parallel = FALSE)
pdf(file="permuTest_outliergenes_0.05_coldspots_less_10000permut.pdf")
plot(Recspt)
dev.off()
lz <- localZScore(A=gr, pt=Recspt, window=10*mean(width(gr)),
                  step=mean(width(gr))/2, B=coldspots)
pdf(file="permuTest_outliergenes_0.05_coldspots_less_10000permut_Zplot.pdf")
plot(lz)
dev.off()

### with masked version of the genome
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10.masked")
genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
genome <- getGenomeAndMask("mm10", mask="BSgenome.Mmusculus.UCSC.mm10.masked")$genome
Recspt <- overlapPermTest(A=gr, B=coldspots, ntimes = 10000, genome=genome, alternative="auto",force.parallel = FALSE)
pdf(file="permuTest_outliergenes_0.05_coldspots_auto_10000permut_masked.pdf")
plot(Recspt)
dev.off()
