library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


peak <- readPeakFile('GSM803365_hg19_wgEncodeHaibTfbsH1hescNrsfV0416102PkRep2.broadpeak')


covplot(peak, weightCol = 'X128.460')

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)


tagHeatmap(tagMatrix, xlim=c(-3000, 3000), xlab = 'Range from TSS', color="red")

options(connectionObserver = NULL) #fixing a bug that did not allow to use the org.Hs.eg.db package
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")


df_peakAnno <- as.data.frame(peakAnno)
write.csv(df_peakAnno, file = 'annotated peaks REST ChIP-seq.csv')



plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAnnoBar(peakAnno)

plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

library(GO.db)

pathway1 <- enrichPathway(df_peakAnno$geneId)


head(pathway1, 2)

library(ReactomePA)
pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
head(pathway1, 2)

gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)

dotplot(pathway2)


library(clusterProfiler)
enrichment <- enricher(df_peakAnno$geneId)


