# New m6A_analysis with diffBind
# Call peaks with MeRIPtools
library("MeRIPtools")
library(ggsci)
setwd("F:/Projects/Sanna/meRIP-seq/")
samplenames <- c( paste0("WB",1:10)) # only on WB
bam_dir <- "F:/Projects/Sanna/meRIP-seq/bams/"
list.files(bam_dir)
MeRIP <- countReads(samplenames = samplenames, gtf = "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf",
                    bamFolder = "F:/Projects/Sanna/meRIP-seq/bams",
                    outputDir = "results",
                    modification = "m6A",
                    fragmentLength = 100,
                    binSize = 50,
                    strandToKeep = "opposite",
                    paired = TRUE,
                    threads = 8
)
summary(MeRIP)
save(MeRIP,file="results/000_rawcounts_WB.rda")

# Perform peak calling for each sample: use the binomial test to account for library coverage
MeRIP <- callPeakBinomial(MeRIP = MeRIP, min_counts = 15, peak_cutoff_fdr = 0.05, peak_cutoff_oddRatio = 1, threads = 7)
# Define joint peak for quantitative analysis.
MeRIP <- reportJointPeak(MeRIPdata = MeRIP, joint_threshold = 2, threads = 8)
head( jointPeak(MeRIP) )
unique.Joint_peak <- MeRIP@jointPeaks[which(!duplicated(paste(MeRIP@jointPeaks$chr,
                                          MeRIP@jointPeaks$start,MeRIP@jointPeaks$end,sep = ":"))),]

# write bed to be used for MOTIF search and DiffBind
write.table(unique.Joint_peak, file = "results/WB_Joint_peak.bed", sep = "\t", row.names = F, col.names = F, quote = F)
# plot peaks distribution over GRanges using Guitar
plotMetaGene(MeRIP@jointPeaks,gtf = "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf")

# Do it again for other samples
# HB (HIV Naive)
samplenames <- c( paste0("HB",1:9)) # only on HB
MeRIP <- countReads(samplenames = samplenames, gtf = "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf",
                    bamFolder = "F:/Projects/Sanna/meRIP-seq/bams",
                    outputDir = "results",
                    modification = "m6A",
                    fragmentLength = 100,
                    binSize = 50,
                    strandToKeep = "opposite",
                    paired = TRUE,
                    threads = 8
)
summary(MeRIP)
save(MeRIP,file="results/000_rawcounts_HB.rda")

# Perform peak calling for each sample: use the binomial test to account for library coverage
MeRIP <- callPeakBinomial(MeRIP = MeRIP, min_counts = 15, peak_cutoff_fdr = 0.05, peak_cutoff_oddRatio = 1, threads = 8)
# Define joint peak for quantitative analysis.
MeRIP <- reportJointPeak(MeRIPdata = MeRIP, joint_threshold = 2, threads = 8)
head( jointPeak(MeRIP) )
unique.Joint_peak <- MeRIP@jointPeaks[which(!duplicated(paste(MeRIP@jointPeaks$chr,
                                                              MeRIP@jointPeaks$start,MeRIP@jointPeaks$end,sep = ":"))),]

# write bed to be used for MOTIF search and DiffBind
write.table(unique.Joint_peak, file = "results/HB_Joint_peak.bed", sep = "\t", row.names = F, col.names = F, quote = F)
# plot peaks distribution over GRanges using Guitar
plotMetaGene(MeRIP@jointPeaks,gtf = "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf")

# WMB (WT + Meth)
samplenames <- c( paste0("WMB",1:10)) # only on WMB
MeRIP <- countReads(samplenames = samplenames, gtf = "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf",
                    bamFolder = "F:/Projects/Sanna/meRIP-seq/bams",
                    outputDir = "results",
                    modification = "m6A",
                    fragmentLength = 100,
                    binSize = 50,
                    strandToKeep = "opposite",
                    paired = TRUE,
                    threads = 8
)
summary(MeRIP)
save(MeRIP,file="results/000_rawcounts_WMB.rda")

# Perform peak calling for each sample: use the binomial test to account for library coverage
MeRIP <- callPeakBinomial(MeRIP = MeRIP, min_counts = 15, peak_cutoff_fdr = 0.05, peak_cutoff_oddRatio = 1, threads = 8)
# Define joint peak for quantitative analysis.
MeRIP <- reportJointPeak(MeRIPdata = MeRIP, joint_threshold = 2, threads = 8)
head( jointPeak(MeRIP) )
unique.Joint_peak <- MeRIP@jointPeaks[which(!duplicated(paste(MeRIP@jointPeaks$chr,
                                                              MeRIP@jointPeaks$start,MeRIP@jointPeaks$end,sep = ":"))),]

# write bed to be used for MOTIF search and DiffBind
write.table(unique.Joint_peak, file = "results/WMB_Joint_peak.bed", sep = "\t", row.names = F, col.names = F, quote = F)
# plot peaks distribution over GRanges using Guitar
plotMetaGene(MeRIP@jointPeaks,gtf = "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf")

# HMB (HIV + Meth)
samplenames <- c( paste0("HMB",1:9)) # only on HMB
MeRIP <- countReads(samplenames = samplenames, gtf = "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf",
                    bamFolder = "F:/Projects/Sanna/meRIP-seq/bams",
                    outputDir = "results",
                    modification = "m6A",
                    fragmentLength = 100,
                    binSize = 50,
                    strandToKeep = "opposite",
                    paired = TRUE,
                    threads = 8
)
summary(MeRIP)
save(MeRIP,file="results/000_rawcounts_HMB.rda")

# Perform peak calling for each sample: use the binomial test to account for library coverage
MeRIP <- callPeakBinomial(MeRIP = MeRIP, min_counts = 15, peak_cutoff_fdr = 0.05, peak_cutoff_oddRatio = 1, threads = 8)
# Define joint peak for quantitative analysis.
MeRIP <- reportJointPeak(MeRIPdata = MeRIP, joint_threshold = 2, threads = 8)
head( jointPeak(MeRIP) )
unique.Joint_peak <- MeRIP@jointPeaks[which(!duplicated(paste(MeRIP@jointPeaks$chr,
                                                              MeRIP@jointPeaks$start,MeRIP@jointPeaks$end,sep = ":"))),]

# write bed to be used for MOTIF search and DiffBind
write.table(unique.Joint_peak, file = "results/HMB_Joint_peak.bed", sep = "\t", row.names = F, col.names = F, quote = F)
# plot peaks distribution over GRanges using Guitar
plotMetaGene(MeRIP@jointPeaks,gtf = "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf")

# # write bed to be used for MOTIF search
# 
# # Exctract peak count to prepare it for differential methylation test
# #MeRIP <- jointPeakCount(MeRIP)
# # name variable names for differential IP test
# #variable(MeRIP) <- data.frame( genotype = c(rep("WT",10),rep("HIV",9)) )
# # take out joint peaks
# Joint_peak <- reportJointPeak(MeRIP, threads = 6)
# unique.Joint_peak <- Joint_peak@jointPeaks[which(!duplicated(paste(Joint_peak@jointPeaks$chr,Joint_peak@jointPeaks$start,Joint_peak@jointPeaks$end,sep = ":"))),]
# # write bed to be used for MOTIF search
# write.table(unique.Joint_peak, file = "results/Naive_Joint_peak.bed", sep = "\t", row.names = F, col.names = F, quote = F)
# plotMetaGene(Joint_peak@jointPeaks,gtf = "G:/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf")