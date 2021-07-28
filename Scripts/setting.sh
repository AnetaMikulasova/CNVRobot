#!/bin/bash



MASTERS="../Masters/"
SH="./sh/"
R="./R/"

#functions
function extract_from_master {
	echo "${line}" | cut -f`awk -F '\t' -v col=$1 'NR==1{for (i=1; i<=NF; i++) if ($i==col) {print i;exit}}' $2`
}



#basic setting
source ${MASTERS}setting_base.sh
QCTABLE=${QC_DIR}"QCmetrics.tsv"



#CN and LOH segmentation setting
SEGMENTATION=${MASTERS}"setting_segmentation.txt"
${DOS2UNIX} ${SEGMENTATION}

SEGMENTATION_SEPARATOR="I"
x=${SEGMENTATION_SEPARATOR}


#masters and qc table
MASTERSEX=${MASTERS}"sex_test_regions.txt"
	${DOS2UNIX} ${MASTERSEX}


#project master
MASTERPROJECTS=${MASTERS}"master_projects.txt"
	${DOS2UNIX} ${MASTERPROJECTS}


#controls master
MASTERCONTROLS=${MASTERS}"master_controls.txt"
${DOS2UNIX} ${MASTERCONTROLS}


#samples master
MASTERSAMPLES=${MASTERS}"master_samples.txt"
${DOS2UNIX} ${MASTERSAMPLES}


#sample and regions of interest master
MASTERREGOFINT=${MASTERS}"master_reg_of_interest_to_plot.txt"
${DOS2UNIX} ${MASTERREGOFINT}


#separators
BEG="\n\n◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎"
END="◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎\n\n"
COL1="\033[1;31m"
COL2="\033[0m"
COL1_WARN="\033[1;33m"
COL1_INFO="\033[1;32m"










