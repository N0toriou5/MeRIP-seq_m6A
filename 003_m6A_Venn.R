# m6A Venn
setwd("F:/Projects/Sanna/meRIP-seq/")
library(DiffBind)
library(GenomicFeatures)
library(VennDiagram)
library(ggforce)
#gtf <- "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf"
gtf <- "F:/genomes/Rat/Rnor6_HIV/Rnor6_hybrid.gtf"
txdb<-makeTxDbFromGFF(gtf,format = "gtf", organism = "Rattus norvegicus",dbxrefTag="gene_name")
# Annotate Diffbind report from genomic ranges to gene symbols ----
library(annotate)
library(ChIPpeakAnno)
annotation<-toGRanges(txdb,feature="gene")
# annotate HIVm vs WTm
load("results/001_DiffBind_HIVm_vs_WTm.rda")
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=0.05)# th is FDR treshold
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
# annotate HIVm vs WTm
load("results/001_DiffBind_HIV_vs_WT.rda")
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=0.05)# th is FDR treshold
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
# Compare with Venn
# A simple single-set diagram
A <- rownames(HIVm_WTm)
B <- rownames(HIV_WT)
overlap <- calculate.overlap(
  x = list(
    "HIVm_vs_WTm" = A,
    "HIV_vs_WT" = B
  )
)
png("plots/003b_Venn_HIVm_vs_WTm.png",w=1500,h=1500,res=300)
draw.pairwise.venn(area1=length(overlap$a1),area2 = length(overlap$a2), cross.area = length(overlap$a3),
                   category = c("HIVm vs. WTm", "HIV vs. WT"), lty = rep("blank", 
                   2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
                   cat.pos = c(0,0), cat.dist = rep(0.025, 2))
dev.off()


# annotate HIVm vs HIV
load("results/001_DiffBind_HIVm_vs_HIV.rda")
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=0.05)# th is FDR treshold
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
HIVm_HIV<-myres[,c(9:12)]
# annotate HIVm vs WTm
load("results/001_DiffBind_WTm_vs_WT.rda")
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=0.05)# th is FDR treshold
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
WTm_WT<-myres[,c(9:12)]
# Compare with Venn
# A simple single-set diagram
A <- rownames(HIVm_HIV)
B <- rownames(WTm_WT)
overlap <- calculate.overlap(
  x = list(
    "HIVm_vs_WTm" = A,
    "HIV_vs_WT" = B
  )
)
png("plots/003b_Venn_HIV_vs_WT-meth.png",w=1500,h=1500,res=300)
draw.pairwise.venn(area1=length(overlap$a1),area2 = length(overlap$a2), cross.area = length(overlap$a3),
                   category = c("HIVm vs. HIV", "WTm vs. WT"), lty = rep("blank", 
                                                                         2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
                   cat.pos = c(0,0), cat.dist = rep(0.025, 2))
dev.off()


# identify exclusively enriched targets in HIV or WT acute
hivm<-overlap$a1
write.table(hivm,file="results/hiv_acute_genes.txt",quote=F,sep="\t",row.names = F,col.names = F)
wtm<-overlap$a2
write.table(wtm,file="results/wt_acute_genes.txt",quote=F,sep="\t",row.names = F,col.names = F)

# enrichR pipeline
# EnrichR GO search

library(enrichR)
dbs <- listEnrichrDbs()
head(dbs)
#dbs <- "KEGG_2019_Human"
#dbs <- "ENCODE_TF_ChIP-seq_2015"
dbs <- c("WikiPathways_2019_Mouse", "KEGG_2019_Mouse")
# 
HIVup <- rownames(HIVm_HIV[HIVm_HIV$Fold > 0,])
HIVdn <- rownames(HIVm_HIV[HIVm_HIV$Fold < 0,])

eup <- enrichr(HIVup, dbs)
edown <- enrichr(HIVdn, dbs)

up <- eup$WikiPathways_2019_Mouse
down <- edown$WikiPathways_2019_Mouse

up$type <- "up"
down$type <- "down"
up <- up [c(1:10),]
up <- up[order(up$Combined.Score), ]  # sort
down <- down [c(1:10),]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
down<-down[order(down$Combined.Score),]
gos <- rbind(down,up)
gos$Term <- gsub("WP(.*)","",gos$Term)
gos$Term <- factor(gos$Term, levels=gos$Term)
png("plots/003b_HIV_GO_Wiki.png",w=2500,h=1500,res=300)
# Diverging Barcharts
ggplot(gos, aes(x=Term, y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=type), width=.5,position="dodge")  +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="lightblue", "up"="#f8766d")) + 
  labs(subtitle="Combined scores from Wiki pathways", 
       title= "Enriched in HIV") + 
  coord_flip() +
  theme_bw()
dev.off()

# Plot also KEGG
up <- eup$KEGG_2019_Mouse
down <- edown$KEGG_2019_Mouse

up$type <- "up"
down$type <- "down"
up <- up [c(1:10),]
up <- up[order(up$Combined.Score), ]  # sort
down <- down [c(1:10),]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
down<-down[order(down$Combined.Score),]
gos <- rbind(down,up)
#gos$Term <- gsub("WP(.*)","",gos$Term)
gos$Term <- factor(gos$Term, levels=gos$Term)
png("plots/003b_HIV_GO_KEGG.png",w=2500,h=1500,res=300)
# Diverging Barcharts
ggplot(gos, aes(x=Term, y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=type), width=.5,position="dodge")  +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="lightblue", "up"="#f8766d")) + 
  labs(subtitle="Combined scores from KEGG pathways", 
       title= "Enriched in HIV") + 
  coord_flip() +
  theme_bw()
dev.off()
