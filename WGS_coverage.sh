#!/bin/bash
#SBATCH --account=fc_insertion1
#SBATCH --partition=savio3
#SBATCH --time=0-6

ml samtools
ml bwa

mkdir stats

bwa mem -t 32 references/chm13v2.0.fa.gz 01-TrimmedData/KCXZ0001D_R1_Trimmed_paired.fastq.gz 01-TrimmedData/KCXZ0001D_R2_Trimmed_paired.fastq.gz > 02-MappedData/KCXZ0001D_T2T.sam

SIZE=$(samtools view -H 02-MappedData/KCXZ0001D_T2T.sam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}')
samtools depth -a 02-MappedData/KCXZ0001D_T2T.sam | awk -v var="$SIZE" '{sum+=$3} END {print sum/var}' > stats/KCXZ0001D_WGScoverage.txt