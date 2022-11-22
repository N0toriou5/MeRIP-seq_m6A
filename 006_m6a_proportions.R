library(Guitar)
setwd("F:/Projects/Sanna/meRIP-seq/results/")
gtf <- "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf"
rn6<-makeTxDbFromGFF(gtf,
                format="gtf",
                dataSource=NA,
                organism="Rattus norvegicus",
                taxonomyId=NA,
                circ_seqs=NULL,
                chrominfo=NULL,
                miRBaseBuild=NA,
                metadata=NULL,
                dbxrefTag="gene_ID")

GRangesListmapToTranscripts
guit_rn6<-makeGuitarTxdb(rn6)
stBedFiles<-"F:/Projects/Sanna/meRIP-seq/results/exomePeak2_WB/Mod.bed"
sites<-blocks(import(stBedFiles))
sitesGranges<-GRangesListmapToTranscripts(sites, mapFilterTranscript = F,transcripts=guit_rn6$tx$tx)
sitesGRangelist <- samplePoints(list(sites),
                             stSampleNum = 5,
                             stAmblguity = 5,
                             pltTxType = c("mrna","ncrna"),
                             stSampleModle = "Equidistance",
                             mapFilterTranscript = FALSE,
                             guit_rn6)
samplePoints(sitesGRangelist,
             stSampleNum = 5,
             stAmblguity = 5,
             pltTxType = c("mrna","ncrna"),
             stSampleModle = "Equidistance",
             mapFilterTranscript = FALSE,
             guit_rn6)
sitesGRangelist$mrna
sitesGRangelist$ncrna$pointsOverlapTx

# All in one image
# HMB vs. HB peaks
stBedFiles<-c("F:/Projects/Sanna/meRIP-seq/results/exomePeak2_HMB/Mod.bed",
              "F:/Projects/Sanna/meRIP-seq/results/exomePeak2_HB/Mod.bed",
              "F:/Projects/Sanna/meRIP-seq/results/exomePeak2_WMB/Mod.bed",
              "F:/Projects/Sanna/meRIP-seq/results/exomePeak2_WB/Mod.bed")
count <- GuitarPlot(
  #txGenomeVer = guitarTxdb,
  txGTF = gtf,
  stBedFiles = stBedFiles,
  miscOutFilePrefix = NA,
  headOrtail = FALSE,
  enableCI = FALSE,
  mapFilterTranscript = TRUE,
  pltTxType = c("ncrna"),
  stGroupName = c("HMB","HB","WMB","WB"))

ncrna<-count$data
boxplot(ncrna$density~ncrna$group,outline=F)

# Use chip profiler
library(ChIPseeker)
peak<-readPeakFile(stBedFiles[1])
peak
peakAnno <- annotatePeak(stBedFiles[4], #tssRegion=c(-3000, 3000),
                         TxDb=rn6,level="transcript", annoDb="org.Rn.eg.db")
plotAnnoPie(peakAnno)

# try again guitar
stGRangeLists = vector("list", length(stBedFiles))
sitesPoints <- list()
for (i in seq_len(length(stBedFiles))) {
  stGRangeLists[[i]] <- blocks(import(stBedFiles[[i]]))
}
for (i in seq_len(length(stGRangeLists))) {
  sitesPoints[[i]] <- samplePoints(stGRangeLists[i],
                                   stSampleNum = 10,
                                   stAmblguity = 5,
                                   pltTxType = c("mrna"),
                                   stSampleModle = "Equidistance",
                                   mapFilterTranscript = FALSE,
                                   guitarTxdb = guit_rn6)
}
table(sitesPoints[[1]]$mrna$pointsOverlapTx)
x<-as.data.frame(sitesPoints[[1]]$mrna)
