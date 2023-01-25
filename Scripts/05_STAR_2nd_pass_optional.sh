#!/bin/bash
#SBATCH --job-name=STAR_2nd_pass
#SBATCH --output=output/%A_%a_STAR_2nd_pass.out
#SBATCH --error=output/%A_%a_STAR_2nd_pass.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --time=40:00:00
#SBATCH --partition=bioinfo
#SBATCH --array=0-79%20

module load STAR/2.5.3a
module load samtools

PROJECT_DIR="/data/bioinfo/data/heart_transcriptomics/EGA_Heesch_Heart/"
cd $PROJECT_DIR

FASTQ_FOLDER=${PROJECT_DIR}/EGA_download/EGAD00001004394

## Map  mRNA-seq samples
if [ ! -d STAR_mRNA_2nd_pass ]; then mkdir STAR_mRNA_2nd_pass; fi

FILES_R1=($FASTQ_FOLDER/PART*/*mR_R1*fastq.gz)
FILES_R2=($FASTQ_FOLDER/PART*/*mR_R2*fastq.gz)
FILE_FULL_PATH_R1=${FILES_R1[$SLURM_ARRAY_TASK_ID]}
FILE_FULL_PATH_R2=${FILES_R2[$SLURM_ARRAY_TASK_ID]}
FILE_NAME=$(echo $FILE_FULL_PATH_R1 | sed "s/_R1_raw.fastq.gz//" | sed "s/.*\///")
NCPU=$(($(nproc)-1))

STAR --runThreadN $NCPU \
--genomeDir /data/public/GENCODE/v38_GRCh38.p13/ \
--readFilesIn $FILE_FULL_PATH_R1 $FILE_FULL_PATH_R2 \
--outFileNamePrefix STAR_mRNA_2nd_pass/${FILE_NAME}_  \
--outSAMtype BAM SortedByCoordinate  \
--quantMode GeneCounts \
--readFilesCommand zcat \
--limitSjdbInsertNsj 2000000 \
--sjdbFileChrStartEnd STAR_mRNA/EGAR00001805112_hs_lv_001_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805113_hs_lv_002_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805114_hs_lv_003_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805115_hs_lv_004_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805116_hs_lv_005_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805117_hs_lv_006_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805118_hs_lv_007_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805119_hs_lv_008_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805120_hs_lv_009_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805121_hs_lv_010_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805122_hs_lv_011_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805123_hs_lv_012_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805124_hs_lv_013_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805125_hs_lv_014_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805126_hs_lv_015_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805127_hs_lv_016_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805128_hs_lv_017_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805129_hs_lv_018_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805130_hs_lv_019_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805131_hs_lv_020_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805132_hs_lv_021_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805133_hs_lv_022_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805134_hs_lv_023_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805135_hs_lv_024_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805136_hs_lv_025_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805137_hs_lv_026_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805138_hs_lv_027_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805139_hs_lv_028_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805140_hs_lv_029_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805141_hs_lv_030_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805142_hs_lv_031_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805143_hs_lv_032_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805144_hs_lv_033_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805145_hs_lv_034_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805146_hs_lv_035_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805147_hs_lv_036_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805148_hs_lv_037_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805149_hs_lv_038_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805150_hs_lv_039_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805151_hs_lv_040_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805152_hs_lv_041_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805153_hs_lv_042_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805154_hs_lv_043_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805155_hs_lv_044_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805156_hs_lv_045_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805157_hs_lv_046_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805158_hs_lv_047_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805159_hs_lv_048_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805160_hs_lv_049_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805161_hs_lv_050_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805162_hs_lv_051_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805163_hs_lv_052_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805164_hs_lv_053_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805165_hs_lv_054_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805166_hs_lv_055_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805167_hs_lv_056_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805168_hs_lv_057_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805169_hs_lv_058_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805170_hs_lv_059_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805171_hs_lv_060_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805172_hs_lv_061_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805173_hs_lv_062_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805174_hs_lv_063_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805175_hs_lv_064_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805176_hs_lv_065_D_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805177_hs_lv_066_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805178_hs_lv_067_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805179_hs_lv_068_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805180_hs_lv_069_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805181_hs_lv_070_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805182_hs_lv_071_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805183_hs_lv_072_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805184_hs_lv_073_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805185_hs_lv_074_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805186_hs_lv_075_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805187_hs_lv_076_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805188_hs_lv_077_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805189_hs_lv_078_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805190_hs_lv_079_C_mR_SJ.out.tab \
STAR_mRNA/EGAR00001805191_hs_lv_081_C_mR_SJ.out.tab
