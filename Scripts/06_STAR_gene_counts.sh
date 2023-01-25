#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH --partition=defq

PROJECT_DIR="/data/bioinfo/data/CRISPR_screen_for_AAV5_infection/STAR_output"
cd $PROJECT_DIR

for i in `seq 1 100`
do
file=`ls ${i}ReadsPerGene.out.tab`
grep "^E" $file | awk -v OFS='\t' '{print $1,$4}' > /data/bioinfo/data/CRISPR_screen_for_AAV5_infection/STAR_gene_counts/${i}_genecounts
done
