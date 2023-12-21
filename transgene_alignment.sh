#!/bin/bash
#SBATCH --account=fc_insertion1
#SBATCH --partition=savio3
#SBATCH --time=0-12

ml samtools
ml bwa
ml java

# deduplicate reads
bash bbmap/clumpify.sh in1=00-RawData/KCXZ0001D_S13_L004_R1_001.fastq.gz \
in2=00-RawData/KCXZ0001D_S13_L004_R2_001.fastq.gz \
out1=00-RawData/KCXZ0001D_R1.dedup.fastq.gz \
out2=00-RawData/KCXZ0001D_R2.dedup.fastq.gz \
dedupe=t optical=f

mkdir 01-TrimmedData

# trim reads
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 00-RawData/KCXZ0001D_R1.dedup.fastq.gz 00-RawData/KCXZ0001D_R2.dedup.fastq.gz \
01-TrimmedData/KCXZ0001D_R1_Trimmed_paired.fastq.gz 01-TrimmedData/KCXZ0001D_R1_Trimmed_unpaired.fastq.gz \
01-TrimmedData/KCXZ0001D_R2_Trimmed_paired.fastq.gz 01-TrimmedData/KCXZ0001D_R2_Trimmed_unpaired.fastq.gz \
SLIDINGWINDOW:6:30 LEADING:30 TRAILING:30 MINLEN:36

mkdir 02-MappedData

# align reads to transgene reference
bwa mem -t 32 references/CBh_template+rDNAflanks.fa 01-TrimmedData/KCXZ0001D_R1_Trimmed_paired.fastq.gz 01-TrimmedData/KCXZ0001D_R2_Trimmed_paired.fastq.gz > 02-MappedData/KCXZ0001D_transgeneflanks-pe.sam

# filter reads, keeping only pairs where one or both mapped to transgene
samtools view -h -F 4 -f 8 02-MappedData/KCXZ0001D_transgeneflanks-pe.sam -o temp.sam
samtools sort temp.sam -o 02-MappedData/KCXZ0001D_out1.sam
samtools view -h -F 8 -f 4 02-MappedData/KCXZ0001D_transgeneflanks-pe.sam -o temp.sam
samtools sort temp.sam -o 02-MappedData/KCXZ0001D_out2.sam
samtools view -h -F 12 02-MappedData/KCXZ0001D_transgeneflanks-pe.sam -o temp.sam
samtools sort temp.sam -o 02-MappedData/KCXZ0001D_out3.sam
samtools merge -f -o 02-MappedData/KCXZ0001D_transgeneflanks_mappedmates.sam 02-MappedData/KCXZ0001D_out1.sam 02-MappedData/KCXZ0001D_out2.sam 02-MappedData/KCXZ0001D_out3.sam
