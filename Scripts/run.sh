#!/bin/bash

CHECKPOINTS="yes"

#extract CNVRobot path
HOME_PATH=`pwd`
RUN_PATH=`echo "$0" | sed -e "s/$(basename "$0")$//g"`
cd ${RUN_PATH}
cd ..
ROBOT_PATH=`pwd`
SETTING=${ROBOT_PATH}"/Scripts/setting.sh"
cd ${HOME_PATH}

#CHECKPOINT - setting.sh exist and is not empty
#============================================================================
if [[ -f ${SETTING} && -s ${SETTING} ]]; then true; else echo "CHECKPOINT ERROR - "${SETTING}" file not recognized, please check its path"; exit; fi
#============================================================================

# #variables that are same for all captures
# source ${SETTING}

##### (0) check master files, databases etc. to prevent errors in setting as much as possible
#####=======================================================================================================
if [[ ${CHECKPOINTS} == "yes" ]]; then ${ROBOT_PATH}"/Scripts/sh/checkpoints.sh" ${SETTING}; fi
if ls ./error 1> /dev/null 2>&1; then
	rm ./error
	exit
fi
#####=======================================================================================================



##### (1) preparation of intervals, controls and gnomad database
#####=======================================================================================================
${ROBOT_PATH}"/Scripts/sh/run_01_prep.sh" ${SETTING}
if ls ./error 1> /dev/null 2>&1; then
	rm ./error
	exit
fi
#####=======================================================================================================



##### (2) samples preparation and analysis
#####=======================================================================================================
${ROBOT_PATH}"/Scripts/sh/run_02_analysis.sh" ${SETTING}
#####=======================================================================================================
