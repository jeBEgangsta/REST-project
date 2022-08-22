library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)


peak <- readPeakFile('GSM1003451_hg19_wgEncodeBroadHistoneH1hescSirt6Pk.broadpeak')


covplot(peak, weightCol = 'X37.647740')

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)


tagHeatmap(tagMatrix, xlim=c(-3000, 3000), xlab = 'Range from TSS', color="red")

options(connectionObserver = NULL) #fixing a bug that did not allow to use the org.Hs.eg.db package
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")


df_peakAnno <- as.data.frame(peakAnno)
write.csv(df_peakAnno, file = 'annotated peaks SIRT6 ChIP-seq.csv')



plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)

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
de <- list(df_peakAnno$geneId)
yy <- enrichKEGG(de[[1]])
kegg_sirt6 <- as.data.frame(yy)

xx <- enrichGO(de[[1]], 'org.Hs.eg.db')
go_sirt6 <- as.data.frame(xx)


dotplot(yy)




peakAnnoList <- lapply('GSM1003451_hg19_wgEncodeBroadHistoneH1hescSirt6Pk.broadpeak', annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")


