# 004_GO-Terms and GSEA on m6A

setwd("F:/Projects/Sanna/meRIP-seq/")
library(DiffBind)
library(GenomicFeatures)
library(VennDiagram)
source("F:/Archive/textplot3.R")
#gtf <- "F:/genomes/Rat/ref/Rattus_norvegicus.Rnor_6.0.101.gtf"
gtf<-"F:/genomes/Rat/Rnor6_HIV/Rnor6_hybrid.gtf"
txdb<-makeTxDbFromGFF(gtf,format = "gtf", organism = "Rattus norvegicus",dbxrefTag="gene_name")
# Annotate Diffbind report from genomic ranges to gene symbols ----
library(annotate)
library(ChIPpeakAnno)
annotation<-toGRanges(txdb,feature="gene")
# annotate HIVm vs WTm
load("results/001_DiffBind_HIVm_vs_WTm.rda")
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=1)# th is FDR treshold
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
HIVm_WTm<-myres[,c(9:12)]
HIVm_WTm$stat<- -log10(HIVm_WTm$FDR)*sign(HIVm_WTm$Fold)

#GSEA
## Load the msigdb database 6.2
library(fgsea)
library(msigdbr)
msigdbr_species()
mdf<-msigdbr(species="Rattus norvegicus",category = "C2") # Retrieve all mouse gene sets
length(mdf$gs_name) # Number of associations 494062
length(length(unique(mdf$gs_name))) # Number of pathways
head(mdf)
mlist<-mdf %>% split(x=.$gene_symbol,f=.$gs_name)
plength<-sapply(mlist,length)
max(plength) # 1734, the biggest pathway #11/12/20

# Calculate enrichment HIVm_WTm

signature<-setNames(HIVm_WTm$stat,rownames(HIVm_WTm))
set.seed(1)
gseas<-fgseaMultilevel(pathways=mlist,stats=signature,eps = 0,minSize=5,maxSize=Inf,nproc=6)
gseas<-gseas[order(gseas$pval),]
save(gseas, file = "results/004_m6A_GSEA_HIVm_WTm.rda")

# annotate HIVm vs HIV
load("results/001_DiffBind_HIVm_vs_HIV.rda")
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=1)# th is FDR treshold
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
HIVm_HIV<-myres[,c(9:12)]
HIVm_HIV$stat<- -log10(HIVm_HIV$FDR)*sign(HIVm_HIV$Fold)

# Calculate enrichment HIVm_HIV

signature<-setNames(HIVm_HIV$stat,rownames(HIVm_HIV))
set.seed(1)
gseas<-fgseaMultilevel(pathways=mlist,stats=signature,eps = 0,minSize=5,maxSize=Inf,nproc=8)
gseas<-gseas[order(gseas$pval),]
save(gseas, file = "results/004_m6A_GSEA_HIVm_HIV.rda")

# annotate WTm vs WT
load("results/001_DiffBind_WTm_vs_WT.rda")
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=1)# th is FDR treshold
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
WTm_WT<-myres[,c(9:12)]
WTm_WT$stat<- -log10(WTm_WT$FDR)*sign(WTm_WT$Fold)

# Calculate enrichment HIVm_HIV

signature<-setNames(WTm_WT$stat,rownames(WTm_WT))
set.seed(1)
gseas<-fgseaMultilevel(pathways=mlist,stats=signature,eps = 0,minSize=5,maxSize=Inf,nproc=8)
gseas<-gseas[order(gseas$pval),]
save(gseas, file = "results/004_m6A_GSEA_WTm_WT.rda")



# Pathway correlation
load("results/004_m6A_GSEA_HIVm_HIV.rda")
hiv<-gseas
load("results/004_m6A_GSEA_WTm_WT.rda")
wt<-gseas
x<-setNames(hiv$NES,hiv$pathway)
y<-setNames(wt$NES,wt$pathway)
common<-intersect(names(x),names(y))
x<-x[common]
y<-y[common]
png("plots/004_scatter_HIVm_vs_HIV.png",w=6500,h=2500,res=300)
plot(x,y,pch=20,col="grey",xlab="HIVm enriched pathways",ylab="WTm enriched pathways",main="Correlation between changes in m6A HIVm vs. WTm at pathway level")
pcc<-cor.test(x,y)
mtext(paste0("R=",signif(pcc$estimate,2)," p=",signif(pcc$p.value,3)))
lml<-lm(y~x)
abline(lml,lwd=1)
#highlight
#x2<-x[grep("REACTOME",names(x),value=F)]
#y2<-y[grep("REACTOME",names(y),value=F)]
topx<-names(sort(x,dec=T)[1:25])
topy<-names(sort(y,dec=T)[1:25])
topup<-intersect(topx,topy)
topx<-names(sort(x,dec=F)[1:25])
topy<-names(sort(y,dec=F)[1:25])
topdn<-intersect(topx,topy)
top<-unique(c(topup,topdn))
#top<-top[-2]
#red<-"SP140"
set.seed(1)
points(x[topup],y[topup],col="red",pch=20)
points(x[topdn],y[topdn],col="blue",pch=20)
textplot3(x[top],y[top],words=top,font=2,cex=0.5,show.lines=F,col="black")
dev.off()

# annotate HIV vs WT
load("results/001_DiffBind_HIV_vs_WT.rda")
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=1)# th is FDR treshold
report<-annotatePeakInBatch(report,featureType = "Exon",
                            AnnotationData=annotation,
                            output="overlapping",
                            FeatureLocForDistance="TSS",
                            bindingRegion=NULL,
                            ignore.strand = FALSE,
                            feature="gene")
# Convert Ensembl names into Gene Symbols
tmp<-ens2eg_Rn(report$feature)
symbols<-eg2sym_Rn(tmp)
report$feature<-symbols
myres<-as.data.frame(report)
myres<-myres[!is.na(myres$feature),]
myres<-myres[!duplicated(myres$feature),]
rownames(myres)<-myres$feature
HIV_WT<-myres[,c(9:12)]
HIV_WT$stat<- -log10(HIV_WT$FDR)*sign(HIV_WT$Fold)

# Calculate enrichment HIVm_WTm

signature<-setNames(HIV_WT$stat,rownames(HIV_WT))
set.seed(1)
gseas<-fgseaMultilevel(pathways=mlist,stats=signature,eps = 0,minSize=5,maxSize=Inf,nproc=8)
gseas<-gseas[order(gseas$pval),]
save(gseas, file = "results/004_m6A_GSEA_HIV_WT.rda")