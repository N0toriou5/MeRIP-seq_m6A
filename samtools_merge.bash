# merge bam files for m6a visualization
cd /mnt/f/Projects/Sanna/meRIP-seq/bams
mkdir merged
samtools merge -@ 8 merged/WB.m6a.bam WB1.m6a.bam WB2.m6a.bam WB3.m6a.bam WB4.m6a.bam WB5.m6a.bam WB6.m6a.bam WB7.m6a.bam WB8.m6a.bam WB9.m6a.bam WB10.m6a.bam
samtools merge -@ 8 merged/WMB.m6a.bam WMB1.m6a.bam WMB2.m6a.bam WMB3.m6a.bam WMB4.m6a.bam WMB5.m6a.bam WMB6.m6a.bam WMB7.m6a.bam WMB8.m6a.bam WMB9.m6a.bam WMB10.m6a.bam
samtools merge -@ 8 merged/HB.m6a.bam HB1.m6a.bam HB2.m6a.bam HB3.m6a.bam HB4.m6a.bam HB5.m6a.bam HB6.m6a.bam HB7.m6a.bam HB8.m6a.bam HB9.m6a.bam
samtools merge -@ 8 merged/HMB.m6a.bam HMB1.m6a.bam HMB2.m6a.bam HMB3.m6a.bam HMB4.m6a.bam HMB5.m6a.bam HMB6.m6a.bam HMB7.m6a.bam HMB8.m6a.bam HMB9.m6a.bam

samtools merge -@ 8 merged/WB.input.bam WB1.input.bam WB2.input.bam WB3.input.bam WB4.input.bam WB5.input.bam WB6.input.bam WB7.input.bam WB8.input.bam WB9.input.bam WB10.input.bam
samtools merge -@ 8 merged/WMB.input.bam WMB1.input.bam WMB2.input.bam WMB3.input.bam WMB4.input.bam WMB5.input.bam WMB6.input.bam WMB7.input.bam WMB8.input.bam WMB9.input.bam WMB10.input.bam
samtools merge -@ 8 merged/HB.input.bam HB1.input.bam HB2.input.bam HB3.input.bam HB4.input.bam HB5.input.bam HB6.input.bam HB7.input.bam HB8.input.bam HB9.input.bam
samtools merge -@ 8 merged/HMB.input.bam HMB1.input.bam HMB2.input.bam HMB3.input.bam HMB4.input.bam HMB5.input.bam HMB6.input.bam HMB7.input.bam HMB8.input.bam HMB9.input.bam

# Script 3
# Bash script
# generate viewable format from bam file
igvtools=/home/notorious/igv/IGV_2.8.11/igvtools
cd /mnt/f/Projects/Sanna/meRIP-seq/
samtools_and_igvtools() {
samtools view -F 516 -q 30 -b ./bams/merged/"$i".bam | samtools sort -o ./bams/merged/"$i"_inter.bam #filter reads
samtools rmdup -s ./bams/merged/"$i"_inter.bam ./bams/merged/"$i"_filtered.bam # remove duplicate reads
$igvtools count -z 5 -w 10 -e 0 ./bams/merged/"$i"_filtered.bam ./tdf/"$i".tdf rn6 #generate TDF
}
export -f samtools_and_igvtools
#execution
for i in HB.input HB.m6a WB.input WB.m6a WMB.input WMB.m6a HMB.input HMB.m6a
do
samtools_and_igvtools $ i
done
PATH=$PATH:/home/notorious/programs/homer/.//bin/