setwd("F:/Projects/Sanna/meRIP-seq/")
library("FastaUtils")
fasta<-"results/diffexomePeak_Hm_Wm/diffHIVmWm.fa"
DNA2RNA(file = fasta)
fasta<-"results/diffexomePeak_Hm_H/diffHIVm.fa"
DNA2RNA(file = fasta)
