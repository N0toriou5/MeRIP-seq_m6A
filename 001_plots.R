#BiocManager::install("Guitar")
# http://127.0.0.1:15151/library/Guitar/doc/Guitar-Overview.pdf
library(Guitar)
setwd("F:/Projects/Sanna/meRIP-seq/results/")
stBedFiles<-c("F:/Projects/Sanna/meRIP-seq/results/exomePeak2_WB/Mod.bed",
              "F:/Projects/Sanna/meRIP-seq/results/exomePeak2_HB/Mod.bed")
#gtf <- "H:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf"
gtf <- "F:/genomes/Rat/ref/Rattus_norvegicus.Rnor_6.0.101.gtf"
#txDB<-makeTxDbFromGFF(gtf,format = "gtf", organism = "Rattus norvegicus")
#guitarTxdb <- makeGuitarTxdb(txdb = txDB, txPrimaryOnly = FALSE)
count <- GuitarPlot(
 # txGenomeVer = guitarTxdb,
                    txGTF = gtf,
                    stBedFiles = stBedFiles,
                    miscOutFilePrefix = NA,
                    headOrtail = TRUE,
                    enableCI = FALSE,
                    mapFilterTranscript = TRUE,
                    pltTxType = c("mrna"),
                    stGroupName = c("WB","HB"))
png("F:/Projects/Sanna/meRIP-seq/plots/001_peak_HIV-WT_mrna.png",w=3500,h=2000,res=300)
plot(count) +
  #theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution on mRNA") + xlab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

count <- GuitarPlot(
  #txGenomeVer = guitarTxdb,
  txGTF = gtf,
  stBedFiles = stBedFiles,
  miscOutFilePrefix = NA,
  headOrtail = T,
  enableCI = FALSE,
  mapFilterTranscript = TRUE,
  pltTxType = c("ncrna"),
  stGroupName = c("WB","HB"))

png("F:/Projects/Sanna/meRIP-seq/plots/001_peak_HIV-WT_ncrna.png",w=3500,h=2000,res=300)
plot(count) +
  #theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution on lncRNA") + xlab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# WMB vs. WB peaks
stBedFiles<-c("F:/Projects/Sanna/meRIP-seq/results/exomePeak2_WMB/Mod.bed",
              "F:/Projects/Sanna/meRIP-seq/results/exomePeak2_HMB/Mod.bed")
count <- GuitarPlot(
  #txGenomeVer = guitarTxdb,
  txGTF = gtf,
  stBedFiles = stBedFiles,
  miscOutFilePrefix = NA,
  headOrtail = TRUE,
  enableCI = FALSE,
  mapFilterTranscript = TRUE,
  pltTxType = c("mrna"),
  stGroupName = c("WMB","HMB"))
png("F:/Projects/Sanna/meRIP-seq/plots/001_peak_HMB-WMB_mrna.png",w=3500,h=2000,res=300)
plot(count) +
  #theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution on mRNA") + xlab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

count <- GuitarPlot(
  #txGenomeVer = guitarTxdb,
  txGTF = gtf,
  stBedFiles = stBedFiles,
  miscOutFilePrefix = NA,
  headOrtail = T,
  enableCI = FALSE,
  mapFilterTranscript = TRUE,
  pltTxType = c("ncrna"),
  stGroupName = c("WMB","HMB"))

png("F:/Projects/Sanna/meRIP-seq/plots/001_peak_HMB-WMB_ncrna.png",w=3500,h=2000,res=300)
plot(count) +
  #theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution on lncRNA") + xlab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# WMB vs. WB peaks
stBedFiles<-c("F:/Projects/Sanna/meRIP-seq/results/exomePeak2_WMB/Mod.bed",
              "F:/Projects/Sanna/meRIP-seq/results/exomePeak2_WB/Mod.bed")
count <- GuitarPlot(
  #txGenomeVer = guitarTxdb,
  txGTF = gtf,
  stBedFiles = stBedFiles,
  miscOutFilePrefix = NA,
  headOrtail = TRUE,
  enableCI = FALSE,
  mapFilterTranscript = TRUE,
  pltTxType = c("mrna"),
  stGroupName = c("WMB","WB"))
png("F:/Projects/Sanna/meRIP-seq/plots/001_peak_WMB-WB_mrna.png",w=3500,h=2000,res=300)
plot(count) +
  #theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution on mRNA") + xlab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

count <- GuitarPlot(
  #txGenomeVer = guitarTxdb,
  txGTF = gtf,
  stBedFiles = stBedFiles,
  miscOutFilePrefix = NA,
  headOrtail = T,
  enableCI = FALSE,
  mapFilterTranscript = TRUE,
  pltTxType = c("ncrna"),
  stGroupName = c("WMB","WB"))

png("F:/Projects/Sanna/meRIP-seq/plots/001_peak_WMB-WB_ncrna.png",w=3500,h=2000,res=300)
plot(count) +
  #theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution on lncRNA") + xlab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# HMB vs. HB peaks
stBedFiles<-c("F:/Projects/Sanna/meRIP-seq/results/exomePeak2_HMB/Mod.bed",
              "F:/Projects/Sanna/meRIP-seq/results/exomePeak2_HB/Mod.bed")
count <- GuitarPlot(
  #txGenomeVer = guitarTxdb,
  txGTF = gtf,
  stBedFiles = stBedFiles,
  miscOutFilePrefix = NA,
  headOrtail = TRUE,
  enableCI = FALSE,
  mapFilterTranscript = TRUE,
  pltTxType = c("mrna"),
  stGroupName = c("HMB","HB"))
png("F:/Projects/Sanna/meRIP-seq/plots/001_peak_HMB-HB_mrna.png",w=3500,h=2000,res=300)
plot(count) +
  #theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution on mRNA") + xlab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

count <- GuitarPlot(
  #txGenomeVer = guitarTxdb,
  txGTF = gtf,
  stBedFiles = stBedFiles,
  miscOutFilePrefix = NA,
  headOrtail = T,
  enableCI = FALSE,
  mapFilterTranscript = TRUE,
  pltTxType = c("ncrna"),
  stGroupName = c("HMB","HB"))

png("F:/Projects/Sanna/meRIP-seq/plots/001_peak_HMB-HB_ncrna.png",w=3500,h=2000,res=300)
plot(count) +
  #theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution on lncRNA") + xlab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

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
  headOrtail = TRUE,
  enableCI = FALSE,
  mapFilterTranscript = TRUE,
  pltTxType = c("mrna"),
  stGroupName = c("HMB","HB","WMB","WB"))
png("F:/Projects/Sanna/meRIP-seq/plots/001_peak_All_mrna.png",w=3500,h=2000,res=300)
plot(count) +
  #theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution on mRNA") + xlab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

count <- GuitarPlot(
  #txGenomeVer = guitarTxdb,
  txGTF = gtf,
  stBedFiles = stBedFiles,
  miscOutFilePrefix = NA,
  headOrtail = T,
  enableCI = FALSE,
  mapFilterTranscript = TRUE,
  pltTxType = c("ncrna"),
  stGroupName = c("HMB","HB","WMB","WB"))

png("F:/Projects/Sanna/meRIP-seq/plots/001_peak_All_ncrna.png",w=3500,h=2000,res=300)
plot(count) +
  #theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("Distribution on lncRNA") + xlab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
