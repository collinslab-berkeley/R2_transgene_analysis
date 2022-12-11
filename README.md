# R2 transgene analysis
#### Bioinformatic analysis package for whole-genome sequencing data following R2-mediated transgene addition

## Overview
The contents of this repository allow analysis of R2 retrotransposon-mediated transgene insertions from whole-genome sequencing data. Corresponding WGS data available on the SRA (`accession TBD`). Bash scripts are tailored for use on Savio, UC Berkeley/LBNL's HPC cluster, with a SLURM scheduling system. Python scripts and jupyter notebooks can be run locally.

To process data:
1. Place raw sequencing data in directory named `00-RawData`
2. Download T2T reference genome and construct alignment indices with `bwa`
3. Run `transgene_alignment.sh`
	- Deduplicate and trim reads
	- Align reads to transgene reference, flanked by 28S rDNA
	- Filter reads, keeping only pairs where one or both mapped to transgene reference
	- NB: this script needs to be run before proceeding with subsequent steps
4. Run `WGS_alignment.sh`
	- Calculates WGS coverage
5. Process and filter reads with `process_transgene_reads.py`
	- Usage: `python process_transgene_reads.py transgene_ref input_sam output_csv`
	- Sample usage: `python process_transgene_reads.py references/220416_template+rDNAflanks.fa ZoAl_TCA5.sam ZoAl_TCA5_reads.csv`
6. Analyze and plot data with `transgene_data_analysis.ipynb`

Code covered under MIT License.

### Dependencies
- Python 3.6+
- bwa v0.7.17
- BBMap v38.97
- Trimmomatic v0.39 (and jdk-17.0.2)
- samtools v1.8

Python packages
- `pip install numpy`
- `pip install pandas`
- `pip install scipy`
- `pip install matplotlib`
- `pip install biopython`
- `pip install pysam`