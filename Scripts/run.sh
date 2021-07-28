#!/bin/bash

#variables that are same for all captures
# source "/Users/aneta/Documents/projects/GATK_CN/Scripts/setting.sh" "/Users/aneta/Documents/projects/GATK_CN/Masters/setting_base.txt"
source ./setting.sh

##### (1) preparation of intervals, controls and gnomad database
#####=======================================================================================================
${SH}run_01_prep.sh
#####=======================================================================================================

##### (2) samples preparation and analysis
#####=======================================================================================================
${SH}run_02_analysis.sh
#####=======================================================================================================
