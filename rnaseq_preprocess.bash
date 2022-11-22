

cd /mnt/f/Projects/Sanna/input_meth
# fastqc
mkdir fastqc
for fname in *.fastq.gz
do
root=`basename $fname .fastq.gz`
echo "Doing $root"
fastqc ${root}.fastq.gz -o fastqc/ -t 7
done

# adapter trimming
mkdir trimmed
for fname in *_R1.fastq.gz
do
root=`basename $fname _R1.fastq.gz`
echo "Doing $root"
fastp -i ${root}_R1.fastq.gz -I ${root}_R2.fastq.gz -o trimmed/${root}_R1_trimmed.fastq.gz -O trimmed/${root}_R2_trimmed.fastq.gz
done

# make hisat2 index
cd trimmed/
mkdir hisat_bams
# Hisat index
index=/mnt/g/Sanna/Rat_Genome/hisat_index/Rat

for fname in *_R1_trimmed.fastq.gz
do
cd /mnt/f/Projects/Sanna/input_meth/trimmed
root=`basename $fname _R1_trimmed.fastq.gz`
echo "Doing $root"
hisat2 -x $index \
-p 7 \
-t \
-1 ${root}_R1_trimmed.fastq.gz \
-2 ${root}_R2_trimmed.fastq.gz \
| samtools view -bS - > hisat_bams/${root}.bam
# Sort and index BAMs
cd /mnt/f/Projects/Sanna/input_meth/trimmed/hisat_bams
samtools sort -m 2G -@ 3 -O BAM -o ${root}.sorted.bam ${root}.bam 
samtools index ${root}.sorted.bam
rm ${root}.bam
done

### Convert to counts with the subread package
cd /mnt/f/Projects/Sanna/meRIP-seq/bams
mkdir rawcounts
gtf=/mnt/g/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf
featureCounts -T 8 -p -s 2 -t exon -g gene_id -a $gtf -o rawcounts/featureCounts.txt *input.bam
featureCounts -T 8 -p -s 2 -t gene -g gene_id -a $gff -o rawcounts/featureCounts_input.txt *input.bam
# call mRNA and ncRNA in separate files
featureCounts=/home/notorious/programs/subread/subread-2.0.1-source/bin/featureCounts
gff=/mnt/g/Sanna/Rat_Genome/GFF/Rattus_norvegicus.Rnor_6.0.100.gff3
featureCounts -T 8 -p -s 2 -t ncRNA_gene -g gene_id -a $gff -o rawcounts/featureCounts_nc.txt *m6A.bam
featureCounts -T 8 -p -s 2 -t gene -g Name -a $gff -o rawcounts/featureCounts.txt *m6A.bam
featureCounts -T 8 -p -s 2 -t exon -g gene_id -a $gff -o rawcounts/featureCounts_exon.txt *m6A.bam

#do it with htseq-count
mkdir counts
for fname in *.m6A.bam
do
root=`basename $fname .m6A.bam`
echo "Doing $root"
htseq-count -s reverse -f bam -t gene -i Name -n 8 ${root}.m6A.bam $gff > counts/${root}.counts