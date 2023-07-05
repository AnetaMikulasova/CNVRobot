#!/bin/bash

CHECKPOINTS="yes"

#extract CNVRobot path
RUN_PATH=`pwd`
SETTING=${RUN_PATH}"/setting.sh"


#CHECKPOINT - setting.sh exist and is not empty
#============================================================================
if [[ -f ${SETTING} && -s ${SETTING} ]]; then true; else echo "CHECKPOINT ERROR - "${SETTING}" file not recognized, please check its path"; exit; fi
#============================================================================



##### (0) check master files, databases etc. to prevent errors in setting as much as possible
#####=======================================================================================================
if [[ ${CHECKPOINTS} == "yes" ]]; then ${RUN_PATH}"/sh/checkpoints.sh" ${SETTING}; fi
if ls ./error 1> /dev/null 2>&1; then
	rm ./error
	exit
fi
#####=======================================================================================================



##### (1) preparation of intervals, controls and gnomad database
#####=======================================================================================================
${RUN_PATH}"/sh/run_01_prep.sh" ${SETTING} ${CHECKPOINTS}
if ls ./error 1> /dev/null 2>&1; then
	rm ./error
	exit
fi
#####=======================================================================================================



##### (2) samples preparation and analysis
#####=======================================================================================================
${RUN_PATH}"/sh/run_02_analysis.sh" ${SETTING} ${CHECKPOINTS}
#####=======================================================================================================
