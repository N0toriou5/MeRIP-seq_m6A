# Differential Binding
# Now analyze for differential m6A_mRNA with DiffBind
# http://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
library(DiffBind)
setwd("F:/Projects/Sanna/meRIP-seq/")
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
save(dba,file="results/rawdba.rda")
png("plots/001_reads_heatmap.png",w=3000,h=3000,res=300)
plot(dba) # Correlation heatmap, using affinity (read count) data
dev.off()

# Classi type PCA
#PCA
png("plots/001_counts_pca.png",w=3000,h=3000,res=300)
dba.plotPCA(dba,label = DBA_ID)
dev.off()

# Establish compare contrasts # do it case x case
# HIVm vs WTm
DBdata <- dba.contrast(dba, group1 = dba$masks$HIVm,group2=dba$masks$WTm,
                       categories=DBA_CONDITION, minMembers = 2,
                       name1 = "HIVm",name2 = "WTm")
DBdata
dba.show(DBdata,bContrasts=T)
# perform differential binding analysis
DBdata<-dba.analyze(DBdata,method=c(DBA_EDGER,DBA_DESEQ2))
DBdata
save(DBdata,file="results/001_DiffBind_HIVm_vs_WTm.rda")
dba.report(DBdata, contrast=1,method=DBA_EDGER)
dba.report(DBdata, contrast=1,method=DBA_DESEQ2)

# HIVm vs HIV
DBdata <- dba.contrast(dba, group1 = dba$masks$HIVm,group2=dba$masks$HIV,
                       categories=DBA_CONDITION, minMembers = 2,
                       name1 = "HIVm",name2 = "HIV")
DBdata
dba.show(DBdata,bContrasts=T)
# perform differential binding analysis
DBdata<-dba.analyze(DBdata,method=c(DBA_EDGER,DBA_DESEQ2))
DBdata
save(DBdata,file="results/001_DiffBind_HIVm_vs_HIV.rda")

# HIV vs WT
DBdata <- dba.contrast(dba, group1 = dba$masks$HIV,group2=dba$masks$WT,
                       categories=DBA_CONDITION, minMembers = 2,
                       name1 = "HIV",name2 = "WT")
DBdata
dba.show(DBdata,bContrasts=T)
# perform differential binding analysis
DBdata<-dba.analyze(DBdata,method=c(DBA_EDGER,DBA_DESEQ2))
DBdata
save(DBdata,file="results/001_DiffBind_HIV_vs_WT.rda")

# WTm vs WT
DBdata <- dba.contrast(dba, group1 = dba$masks$WTm,group2=dba$masks$WT,
                       categories=DBA_CONDITION, minMembers = 2,
                       name1 = "WTm",name2 = "WT")
DBdata
dba.show(DBdata,bContrasts=T)
# perform differential binding analysis
DBdata<-dba.analyze(DBdata,method=c(DBA_EDGER,DBA_DESEQ2))
DBdata
save(DBdata,file="results/001_DiffBind_WTm_vs_WT.rda")
