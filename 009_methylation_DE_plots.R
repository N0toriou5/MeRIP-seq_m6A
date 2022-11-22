library("DiffBind")
load("G:/Daniele/Sanna/meRIP-seq/results/001_DiffBind_HIVm_vs_WTm.rda")
report<-dba.report(DBdata, contrast=1,method=DBA_EDGER,th=1)
library(GenomicFeatures)
gtf <- "F:/genomes/Rat/Rnor6_HIV/Rnor6_hybrid.gtf"
txdb<-makeTxDbFromGFF(gtf,format = "gtf", organism = "Rattus norvegicus",dbxrefTag="gene_name")
# Annotate Diffbind report from genomic ranges to gene symbols ----
library(annotate)
library(ChIPpeakAnno)
annotation<-toGRanges(txdb,feature="gene")
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

# Volcano plot of differential m6A
# Enhanced Volcano
library(EnhancedVolcano)
library(airway)
library("magrittr")            
png("F:/009_m6a_pseudobulk_meth.png",w=2500,h=2500,res=300)
EnhancedVolcano(m6a, subtitle = "",
                lab=rownames(m6a),
                x='Fold',
                y='FDR',
                title="HIVm vs. WTm",
                pCutoff = 0.05,
                FCcutoff = 1,
                labFace = "bold",
                labSize = 4,
                caption = paste0('Enhanced methylation = ', nrow(m6a[m6a$Fold>1&m6a$FDR<=0.05,]), ' genes',"\n",'Reduced methylation = ',
                                 nrow(m6a[m6a$Fold< -0.5&m6a$FDR<=0.05,]), ' genes'),)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()
