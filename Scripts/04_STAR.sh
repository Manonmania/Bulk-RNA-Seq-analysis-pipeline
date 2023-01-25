#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH --partition=highmem

spack load star@2.7.6a
spack load samtools@1.9

for i in `seq 1 192`
do

        fq1=`ls /data/bioinfo/data/CRISPR_screen_for_AAV5_infection/rawREADS/${i}_S*_R1_001.fastq.gz`
        fq2=`ls /data/bioinfo/data/CRISPR_screen_for_AAV5_infection/rawREADS/${i}_S*_R2_001.fastq.gz`

STAR --runThreadN 12 \
--genomeDir /data/bioinfo/data/CRISPR_screen_for_AAV5_infection/STAR_index_EGFP \
--readFilesIn $fq1 $fq2 \
--outFileNamePrefix STAR_output/${i}  \
--outSAMtype BAM SortedByCoordinate  \
--outSAMunmapped Within  \
--outFilterType BySJout  \
--chimSegmentMin 15  \
--chimJunctionOverhangMin 15  \
--chimOutType WithinBAM  \
--twopassMode Basic  \
--quantMode GeneCounts \
--readFilesCommand zcat
done
