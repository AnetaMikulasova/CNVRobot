#!/bin/bash

#extract CNVRobot path
HOME_PATH=`pwd`
RUN_PATH=`echo "$0" | sed -e "s/$(basename "$0")$//g"`
cd ${RUN_PATH}
cd ../..
ROBOT_PATH=`pwd`
cd ${HOME_PATH}

MASTERS=${ROBOT_PATH}"/Masters/"
SH=${ROBOT_PATH}"/Scripts/sh/"
R=${ROBOT_PATH}"/Scripts/R/"

#functions
function extract_from_master {
	echo "${line}" | cut -f`awk -F '\t' -v col=$1 'NR==1{for (i=1; i<=NF; i++) if ($i==col) {print i;exit}}' $2`
}

#separators
BEG="\n\n◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎"
END="◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎◼︎\n\n"
COL1="\033[1;31m"
COL2="\033[0m"
COL1_WARN="\033[1;33m"
COL1_INFO="\033[1;32m"


#CHECKPOINT - master folder exists
#============================================================================
if [[ ${CHECKPOINTS} == "yes" ]]; then
	if [[ -d ${MASTERS} ]]; then true; else 
		echo -e ${COL1}${BEG}${COL2}
		echo -e ${COL1_WARN}"◼︎ CHECKPOINT ERROR - folder not found: "${MASTERS}
		echo -e "◼︎ Please resolve before processing."
		echo -e ${COL1}${END}${COL2}; exit
	fi
fi
#============================================================================


# Master files - define path, check if exist and dos2unix
#============================================================================

#function to check if files exit
function check_master_file_exist {
	if [[ ${CHECKPOINTS} == "yes" ]]; then 
		file=$1
		if [ -f ${file} ]; then true; else
			echo -e ${COL1}${BEG}${COL2}
			echo -e ${COL1_WARN}"◼︎ CHECKPOINT ERROR - master file not found: "${file}
			echo -e "◼︎ Please resolve before processing."
			echo -e ${COL1}${END}${COL2}; exit
		fi
	fi
}

#basic setting
check_master_file_exist ${MASTERS}setting_base.sh
source ${MASTERS}setting_base.sh
QCTABLE=${QC_DIR}"QCmetrics.tsv"

#CN and LOH segmentation setting
SEGMENTATION=${MASTERS}"setting_segmentation.txt"
check_master_file_exist ${SEGMENTATION}
${DOS2UNIX} ${SEGMENTATION}

SEGMENTATION_SEPARATOR="I"
x=${SEGMENTATION_SEPARATOR}

#masters and qc table
MASTERSEX=${MASTERS}"sex_test_regions.txt"
check_master_file_exist ${MASTERSEX}
${DOS2UNIX} ${MASTERSEX}

#project master
MASTERPROJECTS=${MASTERS}"master_projects.txt"
check_master_file_exist ${MASTERPROJECTS}
${DOS2UNIX} ${MASTERPROJECTS}

#controls master
MASTERCONTROLS=${MASTERS}"master_controls.txt"
check_master_file_exist ${MASTERCONTROLS}
${DOS2UNIX} ${MASTERCONTROLS}


#samples master
MASTERSAMPLES=${MASTERS}"master_samples.txt"
check_master_file_exist ${MASTERSAMPLES}
${DOS2UNIX} ${MASTERSAMPLES}

#sample and regions of interest master
MASTERREGOFINT=${MASTERS}"master_reg_of_interest_to_plot.txt"
check_master_file_exist ${MASTERREGOFINT}
${DOS2UNIX} ${MASTERREGOFINT}

#============================================================================











