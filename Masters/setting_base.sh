#!/bin/bash

#1.1 softwares
GATK="/path/to/gatk-folder/gatk"
PICARD="java -jar /path/to/picard-folder/picard.jar"
SAMTOOLS=samtools
BCFTOOLS=bcftools
BEDTOOLS=bedtools
BEDGRAPHTOBIGWIG=bedGraphToBigWig
DOS2UNIX=dos2unix
GUNZIP=gunzip
BGZIP=bgzip
TABIX=tabix

#1.2 folders
RESULTS_DIR=/path/to/CNVRobot/Processing/Results/
SUPPORT_FILES_DIR=/path/to/CNVRobot/Processing/Support_Files/
QC_DIR=/path/to/CNVRobot/Processing/QC_table/
DATABASES_DIR=/path/to/CNVRobot/Databases/

#1.3 noisy SNPs parameters
AFDIF=0.15
AFDEPTH=10
AFPERC=0.05

#NOTES:
# AFDIF = 0.15 — Allele Frequency Difference threshold. The maximum allowed deviation from the expected 0.5 allele fraction (for a true heterozygous SNP) before a SNP is flagged as noisy/unreliable. i.e., only SNPs with AF within 0.5 ± 0.15 (so 0.35–0.65) are kept as "clean" hets.
# AFDEPTH = 10 — Minimum read depth required at a SNP site before its allele fraction is trusted. Below this threshold, AF estimates are noisy due to low sampling, so the SNP gets excluded.
# AFPERC = 0.05 — ??

#1.4 Minimal population allelic frequency of SNP to be selected for allelic count
AF_GNOMAD=0.001


#1.5 CNVkit
#something like "/home/username/anaconda3/etc/profile.d/conda.sh"
CONDA="/path/to/conda.sh"
CNVKIT_ENV="name_of_cnvkit_python_environment"

#1.6 ASCAT
#fraction of SNPs to be included in the analysis
#it should be 1, but might take very long time,
#so for training purpose, good to decrease to some small value like 0.1 
#note that smaller fraction means smaller analysis resolution
ASCAT_FRACTION=1

