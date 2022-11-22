# 002 Compare m6a with RNA-seq
library(DiffBind)
#gtf <- "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf"
gtf <- "F:/genomes/Rat/Rnor6_HIV/Rnor6_hybrid.gtf"
setwd("F:/Projects/Sanna/meRIP-seq/")
source("F:/Archive/textplot3.R")
library(GenomicFeatures)
txdb<-makeTxDbFromGFF(gtf,format = "gtf", organism = "Rattus norvegicus",dbxrefTag="gene_name")
# Annotate Diffbind report from genomic ranges to gene symbols ----
library(annotate)
library(ChIPpeakAnno)
annotation<-toGRanges(txdb,feature="gene")

# annotate HIVm vs WTm
load("results/001_DiffBind_HIVm_vs_WTm.rda")
DBdata
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=1)# th is FDR treshold
# #report<-annotatePeakInBatch(report,
#                             AnnotationData=annotation,
#                             FeatureLocForDistance="middle",
#                             ignore.strand = FALSE)
# #report$feature<-getSYMBOL(report$feature,data='org.Rn.eg')

# try to use precast databases instead of custom txdb
# library(org.Rn.eg.db)
# library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
# annotation<-toGRanges(TxDb.Rnorvegicus.UCSC.rn6.refGene,feature="gene")
report<-annotatePeakInBatch(report,featureType = "Exon",
                            AnnotationData=annotation,
                            output="overlapping",
                            FeatureLocForDistance="TSS",
                            bindingRegion=NULL,
                            ignore.strand = FALSE,
                            feature="gene")
# Convert Ensembl names into Gene Symbols
source("F:/Archive/functions/geneids.R")
tmp<-ens2eg_Rn(report$feature)
symbols<-eg2sym_Rn(tmp)
report$feature<-symbols
myres<-as.data.frame(report)
myres<-myres[!is.na(myres$feature),]
myres<-myres[!duplicated(myres$feature),]
rownames(myres)<-myres$feature
m6a<-myres[,c(9:12)]
load("F:/Projects/Sanna/meRIP-seq/R_Analysis_input/results/001_deres_Meth_HIV_vs_WT.rda")
rna<-results

# scatterplot of differential methylation over input
# plot it
x<-setNames(m6a$Fold,rownames(m6a))
y<-setNames(rna$logFC,rownames(rna))
common<-intersect(names(x),names(y))
x<-x[common]
y<-y[common]
png("plots/002_scatter_HIVm_vs_WTm.png",w=2500,h=2500,res=300)
plot(x,y,pch=20,col="grey",xlab="m6A-mRNA logFC",ylab="mRNA logFC",main="Correlation between changes in HIVm vs. WTm")
pcc<-cor.test(x,y)
mtext(paste0("R=",signif(pcc$estimate,2)," p=",signif(pcc$p.value,3)))
lml<-lm(y~x)
abline(lml,lwd=1)
#highlight
topx<-names(sort(x,dec=T)[1:50])
topy<-names(sort(y,dec=T)[1:50])
topup<-intersect(topx,topy)
topx<-names(sort(x,dec=F)[1:50])
topy<-names(sort(y,dec=F)[1:50])
topdn<-intersect(topx,topy)
top<-unique(c(topup,topdn))
#top<-top[-2]
#red<-"SP140"
set.seed(1)
points(x[topup],y[topup],col="red",pch=20)
points(x[topdn],y[topdn],col="blue",pch=20)
textplot3(x[top],y[top],words=top,font=2,cex=1.2,show.lines=F,col="black")
dev.off()

# Analysis of HIV vs WT contrast
# annotate HIV vs WT
load("results/001_DiffBind_HIV_vs_WT.rda")
DBdata
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=1)# th is FDR treshold
# #report<-annotatePeakInBatch(report,
#                             AnnotationData=annotation,
#                             FeatureLocForDistance="middle",
#                             ignore.strand = FALSE)
# #report$feature<-getSYMBOL(report$feature,data='org.Rn.eg')

# try to use precast databases instead of custom txdb
# library(org.Rn.eg.db)
# library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
# annotation<-toGRanges(TxDb.Rnorvegicus.UCSC.rn6.refGene,feature="gene")
report<-annotatePeakInBatch(report,featureType = "Exon",
                            AnnotationData=annotation,
                            output="overlapping",
                            FeatureLocForDistance="TSS",
                            bindingRegion=NULL,
                            ignore.strand = FALSE,
                            feature="gene")
# Convert Ensembl names into Gene Symbols
source("F:/Archive/functions/geneids.R")
tmp<-ens2eg_Rn(report$feature)
symbols<-eg2sym_Rn(tmp)
report$feature<-symbols
myres<-as.data.frame(report)
myres<-myres[!is.na(myres$feature),]
myres<-myres[!duplicated(myres$feature),]
rownames(myres)<-myres$feature
m6a<-myres[,c(9:12)]
load("F:/Projects/Sanna/meRIP-seq/R_Analysis_input/results/001_deres_Naive_HIV_vs_WT.rda")
rna<-results

# scatterplot of differential methylation over input
# plot it
x<-setNames(m6a$Fold,rownames(m6a))
y<-setNames(rna$logFC,rownames(rna))
common<-intersect(names(x),names(y))
x<-x[common]
y<-y[common]
png("plots/002_scatter_HIV_vs_WT.png",w=2500,h=2500,res=300)
plot(x,y,pch=20,col="grey",xlab="m6A-mRNA logFC",ylab="mRNA logFC",main="Correlation between changes in HIV vs. WT")
pcc<-cor.test(x,y)
mtext(paste0("R=",signif(pcc$estimate,2)," p=",signif(pcc$p.value,3)))
lml<-lm(y~x)
abline(lml,lwd=1)
#highlight
topx<-names(sort(x,dec=T)[1:50])
topy<-names(sort(y,dec=T)[1:50])
topup<-intersect(topx,topy)
topx<-names(sort(x,dec=F)[1:50])
topy<-names(sort(y,dec=F)[1:50])
topdn<-intersect(topx,topy)
top<-unique(c(topup,topdn))
#top<-top[-2]
#red<-"SP140"
set.seed(1)
points(x[topup],y[topup],col="red",pch=20)
points(x[topdn],y[topdn],col="blue",pch=20)
textplot3(x[top],y[top],words=top,font=2,cex=1.2,show.lines=F,col="black")
dev.off()

# Contrast WTm vs WT
# annotate WTm vs WT
load("results/001_DiffBind_WTm_vs_WT.rda")
DBdata
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=1)# th is FDR treshold
# #report<-annotatePeakInBatch(report,
#                             AnnotationData=annotation,
#                             FeatureLocForDistance="middle",
#                             ignore.strand = FALSE)
# #report$feature<-getSYMBOL(report$feature,data='org.Rn.eg')

# try to use precast databases instead of custom txdb
# library(org.Rn.eg.db)
# library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
# annotation<-toGRanges(TxDb.Rnorvegicus.UCSC.rn6.refGene,feature="gene")
report<-annotatePeakInBatch(report,featureType = "Exon",
                            AnnotationData=annotation,
                            output="overlapping",
                            FeatureLocForDistance="TSS",
                            bindingRegion=NULL,
                            ignore.strand = FALSE,
                            feature="gene")
# Convert Ensembl names into Gene Symbols
source("F:/Archive/functions/geneids.R")
tmp<-ens2eg_Rn(report$feature)
symbols<-eg2sym_Rn(tmp)
report$feature<-symbols
myres<-as.data.frame(report)
myres<-myres[!is.na(myres$feature),]
myres<-myres[!duplicated(myres$feature),]
rownames(myres)<-myres$feature
m6a<-myres[,c(9:12)]
load("F:/Projects/Sanna/meRIP-seq/R_Analysis_input/results/001_deres_Meth_WT_vs_WT.rda")
rna<-results

# scatterplot of differential methylation over input
# plot it
x<-setNames(m6a$Fold,rownames(m6a))
y<-setNames(rna$logFC,rownames(rna))
common<-intersect(names(x),names(y))
x<-x[common]
y<-y[common]
png("plots/002_scatter_WTm_vs_WT.png",w=2500,h=2500,res=300)
plot(x,y,pch=20,col="grey",xlab="m6A-mRNA logFC",ylab="mRNA logFC",main="Correlation between changes in WTm vs. WT")
pcc<-cor.test(x,y)
mtext(paste0("R=",signif(pcc$estimate,2)," p=",signif(pcc$p.value,3)))
lml<-lm(y~x)
abline(lml,lwd=1)
#highlight
topx<-names(sort(x,dec=T)[1:50])
topy<-names(sort(y,dec=T)[1:50])
topup<-intersect(topx,topy)
topx<-names(sort(x,dec=F)[1:50])
topy<-names(sort(y,dec=F)[1:50])
topdn<-intersect(topx,topy)
top<-unique(c(topup,topdn))
#top<-top[-2]
#red<-"SP140"
set.seed(1)
points(x[topup],y[topup],col="red",pch=20)
points(x[topdn],y[topdn],col="blue",pch=20)
textplot3(x[top],y[top],words=top,font=2,cex=1.2,show.lines=F,col="black")
dev.off()

# Contrast HIVm vs HIV
# annotate HIVm vs HIV
load("results/001_DiffBind_HIVm_vs_HIV.rda")
DBdata
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=1)# th is FDR treshold
# #report<-annotatePeakInBatch(report,
#                             AnnotationData=annotation,
#                             FeatureLocForDistance="middle",
#                             ignore.strand = FALSE)
# #report$feature<-getSYMBOL(report$feature,data='org.Rn.eg')

# try to use precast databases instead of custom txdb
# library(org.Rn.eg.db)
# library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
# annotation<-toGRanges(TxDb.Rnorvegicus.UCSC.rn6.refGene,feature="gene")
report<-annotatePeakInBatch(report,featureType = "Exon",
                            AnnotationData=annotation,
                            output="overlapping",
                            FeatureLocForDistance="TSS",
                            bindingRegion=NULL,
                            ignore.strand = FALSE,
                            feature="gene")
# Convert Ensembl names into Gene Symbols
source("F:/Archive/functions/geneids.R")
tmp<-ens2eg_Rn(report$feature)
symbols<-eg2sym_Rn(tmp)
report$feature<-symbols
myres<-as.data.frame(report)
myres<-myres[!is.na(myres$feature),]
myres<-myres[!duplicated(myres$feature),]
rownames(myres)<-myres$feature
m6a<-myres[,c(9:12)]
load("F:/Projects/Sanna/meRIP-seq/R_Analysis_input/results/001_deres_Meth_HIV_vs_HIV.rda")
rna<-results

# scatterplot of differential methylation over input
# plot it
x<-setNames(m6a$Fold,rownames(m6a))
y<-setNames(rna$logFC,rownames(rna))
common<-intersect(names(x),names(y))
x<-x[common]
y<-y[common]
png("plots/002_scatter_HIVm_vs_HIV.png",w=2500,h=2500,res=300)
plot(x,y,pch=20,col="grey",xlab="m6A-mRNA logFC",ylab="mRNA logFC",main="Correlation between changes in HIVm vs. HIV")
pcc<-cor.test(x,y)
mtext(paste0("R=",signif(pcc$estimate,2)," p=",signif(pcc$p.value,3)))
lml<-lm(y~x)
abline(lml,lwd=1)
#highlight
topx<-names(sort(x,dec=T)[1:50])
topy<-names(sort(y,dec=T)[1:50])
topup<-intersect(topx,topy)
topx<-names(sort(x,dec=F)[1:50])
topy<-names(sort(y,dec=F)[1:50])
topdn<-intersect(topx,topy)
top<-unique(c(topup,topdn))
#top<-top[-2]
#red<-"SP140"
set.seed(1)
points(x[topup],y[topup],col="red",pch=20)
points(x[topdn],y[topdn],col="blue",pch=20)
textplot3(x[top],y[top],words=top,font=2,cex=1.2,show.lines=F,col="black")
dev.off()
rna["Enpp2",]
