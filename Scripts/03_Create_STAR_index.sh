#!/bin/bash
#SBATCH --job-name=STAR_index
#SBATCH --partition=highmem

module load STAR/2.7.9a 

dir=/data/bioinfo/data/CRISPR_screen_for_AAV5_infection/STAR_index_EGFP
fasta=GRCh38.primary_assembly.genome.fa
gtf=gencode.v38.annotation.gtf

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $dir --genomeFastaFiles $dir/$fasta --sjdbGTFfile $dir/$gtf
