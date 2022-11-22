# Same analysis with exomepeak2 Package
# http://www.bioconductor.org/packages/release/bioc/vignettes/exomePeak2/inst/doc/Vignette_V_0.99.html
library(BiocParallel)
BiocParallel::register(BiocParallel::SnowParam(8))
library("exomePeak2")
setwd("F:/Projects/Sanna/meRIP-seq/")
gtf <- "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf"
#genome <- "F:/genomes/Rat/FASTA/rn6.fa"
m6A_filenames<-dir(path="F:/Projects/Sanna/meRIP-seq/bams",pattern="*.m6A.bam$",full.names = T)
input_filenames<-dir("F:/Projects/Sanna/meRIP-seq/bams",pattern=".input.bam$",full.names = T)
m6A<-grep("WB",m6A_filenames,value=T)
input<-grep("WB",input_filenames,value=T)
##### Single samples analysis for WB
peaks<-exomePeak2(bam_ip = m6A,
           bam_input = input,
           gff_dir = gtf,
           paired_end = TRUE,
           library_type = "2nd_strand",
           parallel = TRUE,
           fragment_length = 100,
           save_dir = "exomePeak2_WB",
           peak_calling_mode = "exon")
save(peaks,file="results/exome2_peaks_WB.rda")

# HB (Naive HIV)
m6A<-grep("HB",m6A_filenames,value=T)
input<-grep("HB",input_filenames,value=T)
##### Single samples analysis for WB
peaks<-exomePeak2(bam_ip = m6A,
                  bam_input = input,
                  gff_dir = gtf,
                  paired_end = TRUE,
                  library_type = "2nd_strand",
                  parallel = TRUE,
                  fragment_length = 100,
                  save_dir = "exomePeak2_HB",
                  peak_calling_mode = "exon")
save(peaks,file="results/exome2_peaks_HB.rda")

# WMB (WT + Meth)
m6A<-grep("WMB",m6A_filenames,value=T)
input<-grep("WMB",input_filenames,value=T)
##### Single samples analysis for WB
peaks<-exomePeak2(bam_ip = m6A,
                  bam_input = input,
                  gff_dir = gtf,
                  paired_end = TRUE,
                  library_type = "2nd_strand",
                  parallel = TRUE,
                  fragment_length = 100,
                  save_dir = "exomePeak2_WMB",
                  peak_calling_mode = "exon")
save(peaks,file="results/exome2_peaks_WMB.rda")

# HMB (HIV + Meth)
m6A<-grep("HMB",m6A_filenames,value=T)
input<-grep("HMB",input_filenames,value=T)
##### Single samples analysis for WB
peaks<-exomePeak2(bam_ip = m6A,
                  bam_input = input,
                  gff_dir = gtf,
                  paired_end = TRUE,
                  library_type = "2nd_strand",
                  parallel = TRUE,
                  fragment_length = 100,
                  save_dir = "exomePeak2_HMB",
                  peak_calling_mode = "exon")
save(peaks,file="results/exome2_peaks_HMB.rda")

# Now analyze for differential m6A_mRNA with DiffBind
# http://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
library(DiffBind)
samples<-read.csv("results/Samples_DiffBind_2.csv",header=T)
names(samples)
samples

# read in samples to create a DBA object
dba <- dba(sampleSheet=samples)
# cross-correlations clustering plot of each row of the binding matrix:
plot(dba) # Correlation heatmap, using occupancy (peak caller score) data
# count the reads
dba <- dba.count(dba)
