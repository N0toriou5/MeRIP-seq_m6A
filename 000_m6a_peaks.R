# 000 Sanna m6A dataset
library(BiocParallel)
BiocParallel::register(BiocParallel::SnowParam(8))
library("exomePeak2")
setwd("F:/Projects/Sanna/meRIP-seq/")
#gtf <- "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf"
gtf <- "F:/genomes/Rat/Rnor6_HIV/Rnor6_hybrid.gtf"
m6A_filenames<-dir(path="F:/Projects/Sanna/meRIP-seq/bams",pattern="*.m6a.bam$",full.names = T)
input_filenames<-dir("F:/Projects/Sanna/meRIP-seq/bams",pattern=".input.bam$",full.names = T)
#bam_dir <- "F:/Projects/Sanna/meRIP-seq/bams/"
samplenames<-gsub("F:/Projects/Sanna/meRIP-seq/bams/","",m6A_filenames)
samplenames<-gsub(".m6a.bam","",samplenames)
for (i in 1:length(m6A_filenames)){
m6A<-m6A_filenames[i]
input<-input_filenames[i]
##### Single samples analysis for WB
peaks<-exomePeak2(bam_ip = m6A,
                  bam_input = input,
                  gff_dir = gtf,
                  paired_end = TRUE,
                  library_type = "2nd_strand",
                  parallel = TRUE,
                  fragment_length = 100,
                  save_dir = paste0("results/",samplenames[i],"_singlepeaks"),
                  peak_calling_mode = "exon")
save(peaks,file=paste0("results/",samplenames[i],"_singlepeaks.rda"))
}

# Infinite values are produced whenever Input has 0 counts, so IP/Input=Inf
# To deal with it, we put a maximum read count to Inf values, eg.mybed[is.infinite(mybed[,5]),5]<-1000
for (name in samplenames){
  mybed<-read.delim(paste0("results/",name,"_singlepeaks/Mod.bed"),header=F)
  mybed[is.infinite(mybed[,5]),5]<-1000
  write.table(mybed,file=paste0("results/",name,"_singlepeaks/",name,".bed"),col.names=F,row.names = F,quote=F,sep="\t")
}

# Now analyze for differential m6A_mRNA with DiffBind
# http://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
library(DiffBind)
samples<-read.csv("results/Samples_DiffBind_2.csv",header=T)
names(samples)
samples

# read in samples to create a DBA object
dba <- dba(sampleSheet=samples)
# cross-correlations clustering plot of each row of the binding matrix:
png("plots/000_peaks_heatmap.png",w=3000,h=3000,res=300)
plot(dba) # Correlation heatmap, using occupancy (peak caller score) data
dev.off()
#PCA
png("plots/001_peaks_pca.png",w=3000,h=3000,res=300)
dba.plotPCA(dba,label=DBA_ID)
dev.off()

# count the reads
dba <- dba.count(dba)
dba
info <- dba.show(dba)
png("plots/000_reads_heatmap.png",w=3000,h=3000,res=300)
plot(dba) # Correlation heatmap, using affinity (read count) data
dev.off()

# Classi type PCA
#PCA
png("plots/001_counts_pca.png",w=3000,h=3000,res=300)
dba.plotPCA(dba,label = DBA_ID)
dev.off()

# Establish compare contrasts
DBdata <- dba.contrast(dba, group1 = dba$masks$HIVm,group2=dba$masks$WTm,
                       categories=DBA_CONDITION, minMembers = 2,
                       name1 = "HIVm",name2 = "WTm")
DBdata
dba.show(DBdata,bContrasts=T)
# perform differential binding analysis
DBdata<-dba.analyze(DBdata,method=c(DBA_EDGER,DBA_DESEQ2))
DBdata
# normalize the dataset by library size
dba.norm <- dba.normalize(dba)
par(mfrow=c(2,1))
dba.plotMA(DBdata, contrast=1)
dba.plotMA(DBdata,contrast=1, bNormalized=FALSE)                  
dba.report(DBdata, contrast=1,method=DBA_EDGER)
dba.report(DBdata, contrast=1,method=DBA_DESEQ2)
