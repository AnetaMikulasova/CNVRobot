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
RESULTS_DIR=/path/to/CNVRobot/Procesing/Results/
SUPPORT_FILES_DIR=/path/to/CNVRobot/Procesing/Support_Files/
QC_DIR=/path/to/CNVRobot/Procesing/QC_table/
DATABASES_DIR=/path/to/CNVRobot/Databases/

#1.3 noisy SNPs parameters
AFDIF=0.15
AFDEPTH=10
# AFDEPTH=5
AFPERC=0.05

#1.4 minimal allelic frequency of SNP to be selected
AF_GNOMAD=0.001