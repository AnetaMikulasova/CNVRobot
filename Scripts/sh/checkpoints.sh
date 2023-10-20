#!/bin/bash

# source ./setting.sh
source ${1}

echo -e ${COL1}${BEG}${COL2}
echo -e ${COL1}"◼︎ $(date)"${COL2}
echo -e ${COL1}"◼︎ Testing setting ..."${COL2}
echo -e ${COL1}${END}${COL2}


#CHECKPOINT - databases folder and all expected files exist 
#============================================================================
if [[ -d ${DATABASES_DIR} ]]; then true; else
	echo -e ${COL1}${BEG}${COL2}
	echo -e ${COL1_WARN}"◼︎ CHECKPOINT ERROR - database folder not found: "${DATABASES_DIR}
	echo -e "◼︎ Databases are part of the CNVRobot package. See GitHub website for link to download."
	echo -e "◼︎ Please resolve before processing."
	echo -e ${COL1}${END}${COL2}; 
	touch ./error
	exit
fi

ERROR_DTB=0

#function to check if files exit
function check_database_file_exist {
	file=$1
	if [[ -f ${file} ]]; then true; else
		echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - file in databases not found: "${file}
		echo -e ""${COL2}
		((ERROR_DTB=ERROR_DTB+1))
	fi
}

for FILE in \
${DATABASES_DIR}DGV/DGV_CHM13v2.0_liftover_processed_variants_2020-02-25.rds \
${DATABASES_DIR}DGV/DGV_GRCh38hg38_processed_variants_2020-02-25.rds \
${DATABASES_DIR}DGV/DGV_hg19_processed_variants_2020-02-25.rds \
\
${DATABASES_DIR}ENCODE_BLACKLIST/ENCODE_blacklist_chm13v2.0_liftover.txt \
${DATABASES_DIR}ENCODE_BLACKLIST/ENCODE_blacklist_GRCh38_hg38_liftover.txt \
${DATABASES_DIR}ENCODE_BLACKLIST/ENCODE_blacklist_hg19.txt \
\
${DATABASES_DIR}GENE_ANNOTATION/RefSeqAll_hg19_191206.txt \
${DATABASES_DIR}GENE_ANNOTATION/RefSeqAll_hg38_191206.txt \
\
${DATABASES_DIR}GNOMAD/gnomad_2.1_GRCh37/gnomad_2.1_combine_processed_GRCh37_22XY.vcf.gz \
${DATABASES_DIR}GNOMAD/gnomad_2.1_GRCh37/gnomad_2.1_combine_processed_GRCh37_22XY.vcf.gz.tbi \
${DATABASES_DIR}GNOMAD/gnomad_2.1_GRCh38_liftover/gnomad_2.1_combine_processed_GRCh38_22XY.vcf.gz \
${DATABASES_DIR}GNOMAD/gnomad_2.1_GRCh38_liftover/gnomad_2.1_combine_processed_GRCh38_22XY.vcf.gz.tbi \
${DATABASES_DIR}GNOMAD/gnomad_2.1_hg19_liftover/gnomad_2.1_combine_processed_hg19_22XY.vcf.gz \
${DATABASES_DIR}GNOMAD/gnomad_2.1_hg19_liftover/gnomad_2.1_combine_processed_hg19_22XY.vcf.gz.tbi \
${DATABASES_DIR}GNOMAD/gnomad_2.1_hg38_liftover/gnomad_2.1_combine_processed_hg38_22XY.vcf.gz \
${DATABASES_DIR}GNOMAD/gnomad_2.1_hg38_liftover/gnomad_2.1_combine_processed_hg38_22XY.vcf.gz.tbi \
\
${DATABASES_DIR}MAPPABILITY/mappability_hg19.rds \
${DATABASES_DIR}MAPPABILITY/mappability_hg38GRCh38.rds \
\
${DATABASES_DIR}T2T/chm13v2.0/catLiftOffGenesV1_edit.bed \
${DATABASES_DIR}T2T/chm13v2.0/chm13v2.0_100mer_mappability.rds \
${DATABASES_DIR}T2T/chm13v2.0/chm13v2.0_centromeres.txt \
${DATABASES_DIR}T2T/chm13v2.0/chm13v2.0_cytoBand.txt \
${DATABASES_DIR}T2T/chm13v2.0/chm13v2.0_dbSNPv155_gnomad_processed.vcf.gz \
${DATABASES_DIR}T2T/chm13v2.0/chm13v2.0_dbSNPv155_gnomad_processed.vcf.gz.tbi \
\
${DATABASES_DIR}UCSC/chr_bands/GRCh37_cytoBand.txt \
${DATABASES_DIR}UCSC/chr_bands/GRCh38_cytoBand.txt \
${DATABASES_DIR}UCSC/chr_bands/GRCh38_hg38_cytoBand.txt \
${DATABASES_DIR}UCSC/chr_bands/hg19_cytoBand.txt \
${DATABASES_DIR}UCSC/chr_centromere/GRCh37_centromeres.txt \
${DATABASES_DIR}UCSC/chr_centromere/GRCh38_centromeres.txt \
${DATABASES_DIR}UCSC/chr_centromere/GRCh38_hg38_centromeres.txt \
${DATABASES_DIR}UCSC/chr_centromere/hg19_centromeres.txt
do
	check_database_file_exist ${FILE}
done

if [[ ${ERROR_DTB} == 0 ]]; then true; else
	echo -e ${COL1}${BEG}${COL2}
	echo -e ${COL1}"◼︎ $(date)"${COL2}
	echo -e ${COL1_WARN}"◼︎ CHECKPOINT ERROR - total number of "${ERROR_DTB}" files not found in databases folder (listed above)"
	echo -e "◼︎ These files are part of the CNVRobot package. See GitHub website for link to download."
	echo -e "◼︎ Please resolve before processing."
	echo -e ${COL1}${END}${COL2}
	touch ./error
	exit
fi
#============================================================================



#CHECKPOINT - projects master
#============================================================================
#checkpoint for Databases folder exist + all master files exist - in setting.sh

#Is any project selected?
NofPROJECTs_count=`cat ${MASTERPROJECTS} | awk -F"\t" '$1=="yes"' | wc -l`
NofPROJECTs=$(echo $NofPROJECTs_count | tr -d ' ')
if [[ $NofPROJECTs == 0 ]]; then 
	echo -e ${COL1}${BEG}${COL2}
	echo -e ${COL1_WARN}"◼︎ CHECKPOINT ERROR - no project selected in "${MASTERPROJECTS}" to be run"${COL2}
	echo -e ${COL1}${END}${COL2}
	touch ./error
	exit
fi

ERROR_MASTERPROJECT="./error_master"
touch $ERROR_MASTERPROJECT

cat ${MASTERPROJECTS} | awk -F"\t" '$1=="yes"' | while read -r line || [[ -n "$line" ]]; do

	PROJECT_ID=`extract_from_master 'PROJECT_ID' ${MASTERPROJECTS}`
	PROJECT_TYPE=`extract_from_master 'PROJECT_TYPE' ${MASTERPROJECTS}`
	SEQ_TYPE=`extract_from_master 'SEQ_TYPE' ${MASTERPROJECTS}`
	CAPTURE_ID=`extract_from_master 'CAPTURE_ID' ${MASTERPROJECTS}`
	CAPTURE_FILE=`extract_from_master 'CAPTURE_FILE' ${MASTERPROJECTS}`
	GENOME_VERSION=`extract_from_master 'GENOME_VERSION' ${MASTERPROJECTS}`
	REF=`extract_from_master 'REF' ${MASTERPROJECTS}`
	SETTING_MODE=`extract_from_master 'SETTING_MODE' ${MASTERPROJECTS}`
	GENE_ANNOTATION=`extract_from_master 'GENE_ANNOTATION' ${MASTERPROJECTS}`
	BIN=`extract_from_master 'BIN' ${MASTERPROJECTS}`
	PADDING=`extract_from_master 'PADDING' ${MASTERPROJECTS}`
	WAY_TO_BAM=`extract_from_master 'WAY_TO_BAM' ${MASTERPROJECTS}`
	CTRL=`extract_from_master 'CTRL' ${MASTERPROJECTS}`
	CTRL_BAM_DIR=`extract_from_master 'CTRL_BAM_DIR' ${MASTERPROJECTS}`
	CTRL_BAM_PATTERN=`extract_from_master 'CTRL_BAM_PATTERN' ${MASTERPROJECTS}`
	CTRL_SEX_TEST=`extract_from_master 'CTRL_SEX_TEST' ${MASTERPROJECTS}`
	CTRL_PON_SEX_SELECT=`extract_from_master 'CTRL_PON_SEX_SELECT' ${MASTERPROJECTS}`
	SMPL_BAM_DIR=`extract_from_master 'SMPL_BAM_DIR' ${MASTERPROJECTS}`
	SMPL_BAM_PATTERN=`extract_from_master 'SMPL_BAM_PATTERN' ${MASTERPROJECTS}`
	SMPL_PON_SEX_SELECT=`extract_from_master 'SMPL_PON_SEX_SELECT' ${MASTERPROJECTS}`
	CN_FREQ_IN_CTRLS=`extract_from_master 'CN_FREQ_IN_CTRLS' ${MASTERPROJECTS}`
	GNOMAD_SELECTION=`extract_from_master 'GNOMAD_SELECTION' ${MASTERPROJECTS}`
	SEGMENTATION_ID=`extract_from_master 'SEGMENTATION_ID' ${MASTERPROJECTS}`
	SEGMENTATION_ID_USE=`extract_from_master 'SEGMENTATION_ID_USE' ${MASTERPROJECTS}`

	[ ${SETTING_MODE} == "default" ] && SETTING_MODE="F"
	READ_COUNT_FILTER="${SETTING_MODE:0:1}"
	[[ ${READ_COUNT_FILTER} == "-" ]] && READ_COUNT_FILTER="F"



	#master project - check colums does not contain anything unexpected
	#function to check if provided value does not contain space character
	function check_space_in_master {
		master_file=$1
		master_column=$2
		given_value=$3
		project_id=$4
		if [[ "$given_value" = *" "* ]]; then
			echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - space character detected in "${master_file}
			echo -e "   Project: "${project_id}
			echo -e "   Column: "${master_column}
			echo -e "   Given value (in brackets to see all space characters): ("${given_value}")"
			echo -e ""${COL2}
			echo "error" >> $ERROR_MASTERPROJECT
		fi
	}

	check_space_in_master "${MASTERPROJECTS}" "PROJECT_ID" "${PROJECT_ID}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "PROJECT_TYPE" "${PROJECT_TYPE}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "SEQ_TYPE" "${SEQ_TYPE}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "CAPTURE_ID" "${CAPTURE_ID}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "CAPTURE_FILE" "${CAPTURE_FILE}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "GENOME_VERSION" "${GENOME_VERSION}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "REF" "${REF}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "SETTING_MODE" "${SETTING_MODE}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "GENE_ANNOTATION" "${GENE_ANNOTATION}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "BIN" "${BIN}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "PADDING" "${PADDING}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "WAY_TO_BAM" "${WAY_TO_BAM}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "CTRL" "${CTRL}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "CTRL_BAM_DIR" "${CTRL_BAM_DIR}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "CTRL_BAM_PATTERN" "${CTRL_BAM_PATTERN}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "CTRL_SEX_TEST" "${CTRL_SEX_TEST}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "CTRL_PON_SEX_SELECT" "${CTRL_PON_SEX_SELECT}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "SMPL_BAM_DIR" "${SMPL_BAM_DIR}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "SMPL_BAM_PATTERN" "${SMPL_BAM_PATTERN}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "SMPL_PON_SEX_SELECT" "${SMPL_PON_SEX_SELECT}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "CN_FREQ_IN_CTRLS" "${CN_FREQ_IN_CTRLS}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "GNOMAD_SELECTION" "${GNOMAD_SELECTION}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "SEGMENTATION_ID" "${SEGMENTATION_ID}" "${PROJECT_ID}"
	check_space_in_master "${MASTERPROJECTS}" "SEGMENTATION_ID_USE" "${SEGMENTATION_ID_USE}" "${PROJECT_ID}"


	#master project - check colums have anything but shouldn't be empty
	#function to check if provided value is not empty
	function check_empty_in_master {
		master_file=$1
		master_column=$2
		given_value=$3
		project_id=$4
		if [[ -z "$given_value" ]]; then
			echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - unexpected empty value in "${master_file}
			echo -e "   Project: "${project_id}
			echo -e "   Column: "${master_column}
			echo -e ""${COL2}
			echo "error" >> $ERROR_MASTERPROJECT
		fi
	}
	
	check_empty_in_master "${MASTERPROJECTS}" "PROJECT_ID" "${PROJECT_ID}" "${PROJECT_ID}"
	check_empty_in_master "${MASTERPROJECTS}" "CAPTURE_ID" "${CAPTURE_ID}" "${PROJECT_ID}"
	if [[ ${WAY_TO_BAM} == "find_in_dir" ]]; then check_empty_in_master "${MASTERPROJECTS}" "CTRL_BAM_PATTERN" "${CTRL_BAM_PATTERN}" "${PROJECT_ID}"; fi
	if [[ ${WAY_TO_BAM} == "find_in_dir" ]]; then check_empty_in_master "${MASTERPROJECTS}" "SMPL_BAM_PATTERN" "${SMPL_BAM_PATTERN}" "${PROJECT_ID}"; fi
	check_empty_in_master "${MASTERPROJECTS}" "SEGMENTATION_ID" "${SEGMENTATION_ID}" "${PROJECT_ID}"



	#master project - check columns with specific values expected
	#function to check if provided value match expected values
	function check_value_in_master {
		master_file=$1
		master_column=$2
		given_value=$3
		possible_values=$4
		project_id=$5
		if [[ ${possible_values} =~ (^|[[:space:]])"${given_value}"($|[[:space:]]) ]]; then true; else
				echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - unexpected value in "${master_file}
				echo -e "   Project: "${project_id}
				echo -e "   Column: "${master_column}
				echo -e "   Given value: "${given_value}
				echo -e "   Possible values: "${possible_values}" (see README.txt for details)"
				echo -e ""${COL2}
				echo "error" >> $ERROR_MASTERPROJECT
		fi
	}

	PROJECT_TYPE_VALUES=("germline" "tumor" "other")
	SEQ_TYPE_VALUES=("WGS" "WES" "TS")
	GENOME_VERSION_VALUES=("GRCh37-hg19" "GRCh38-hg38" "CHM13v2.0")
	READ_COUNT_FILTER_VALUES=("default" "F" "D")
	WAY_TO_BAM_VALUES=("default" "absolute" "find_in_dir")
	CTRL_VALUES=("yes" "no")
	CTRL_SEX_TEST_VALUES=("default" "custom" "no")
	CTRL_PON_SEX_SELECT_VALUES=("default" "M" "F" "matched" "mixed")
	SMPL_PON_SEX_SELECT_VALUES=("default" "M" "F" "mixed" "both_separated" "matched_main" "matched_each" "all_options")
	CN_FREQ_IN_CTRLS_VALUES=("yes" "no")
	GNOMAD_SELECTION_VALUES=("default" "capture_filter" "capture_filter_and_ctrl_denois" "af_filter" "af_filter_and_ctrl_denois")
	SEGMENTATION_ID_USE_VALUES=("default" "segm_id" "full")
	check_value_in_master "${MASTERPROJECTS}" "PROJECT_TYPE" "${PROJECT_TYPE}" "${PROJECT_TYPE_VALUES[*]}" "${PROJECT_ID}" 
	check_value_in_master "${MASTERPROJECTS}" "SEQ_TYPE" "${SEQ_TYPE}" "${SEQ_TYPE_VALUES[*]}" "${PROJECT_ID}" 
	check_value_in_master "${MASTERPROJECTS}" "GENOME_VERSION" "${GENOME_VERSION}" "${GENOME_VERSION_VALUES[*]}" "${PROJECT_ID}" 
	check_value_in_master "${MASTERPROJECTS}" "SETTING_MODE" "${READ_COUNT_FILTER}" "${READ_COUNT_FILTER_VALUES[*]}" "${PROJECT_ID}" 
	check_value_in_master "${MASTERPROJECTS}" "WAY_TO_BAM" "${WAY_TO_BAM}" "${WAY_TO_BAM_VALUES[*]}" ${PROJECT_ID}
	check_value_in_master "${MASTERPROJECTS}" "CTRL" "${CTRL}" "${CTRL_VALUES[*]}" "${PROJECT_ID}" 
	check_value_in_master "${MASTERPROJECTS}" "CTRL_SEX_TEST" "${CTRL_SEX_TEST}" "${CTRL_SEX_TEST_VALUES[*]}" "${PROJECT_ID}" 
	check_value_in_master "${MASTERPROJECTS}" "CTRL_PON_SEX_SELECT" "${CTRL_PON_SEX_SELECT}" "${CTRL_PON_SEX_SELECT_VALUES[*]}" "${PROJECT_ID}" 
	check_value_in_master "${MASTERPROJECTS}" "SMPL_PON_SEX_SELECT" "${SMPL_PON_SEX_SELECT}" "${SMPL_PON_SEX_SELECT_VALUES[*]}" "${PROJECT_ID}" 
	check_value_in_master "${MASTERPROJECTS}" "CN_FREQ_IN_CTRLS" "${CN_FREQ_IN_CTRLS}" "${CN_FREQ_IN_CTRLS_VALUES[*]}" "${PROJECT_ID}" 
	check_value_in_master "${MASTERPROJECTS}" "GNOMAD_SELECTION" "${GNOMAD_SELECTION}" "${GNOMAD_SELECTION_VALUES[*]}" "${PROJECT_ID}" 
	check_value_in_master "${MASTERPROJECTS}" "SEGMENTATION_ID_USE" "${SEGMENTATION_ID_USE}" "${SEGMENTATION_ID_USE_VALUES[*]}" "${PROJECT_ID}" 


	#master project - check file that should exit
	#function to check if files given in master project exist
	function check_file_exist_in_master {
		master_file=$1
		master_column=$2
		given_file=$3
		project_id=$4
		if [[ -f ${given_file} ]]; then true; else
			echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - could not find file filled in "${master_file}
			echo -e "   Project: "${project_id}
			echo -e "   Column: "${master_column}
			echo -e "   Given file path: "${given_file}
			echo -e ""${COL2}
			echo "error" >> $ERROR_MASTERPROJECT
		fi
	}
	
	if [[ "${SEQ_TYPE}" == "WES" ]] || [[ "${SEQ_TYPE}" == "TS" ]]; then check_file_exist_in_master "${MASTERPROJECTS}" "CAPTURE_FILE" "${CAPTURE_FILE}" "${PROJECT_ID}"; fi
	check_file_exist_in_master "${MASTERPROJECTS}" "REF" "${REF}" "${PROJECT_ID}"
	if [[ ${GENE_ANNOTATION} != "default" ]]; then check_file_exist_in_master "${MASTERPROJECTS}" "GENE_ANNOTATION" "${GENE_ANNOTATION}" "${PROJECT_ID}"; fi

	#master project - check values expected to be numerical
	#function to check if values values are numerical
	function check_value_num_in_master {
		master_file=$1
		master_column=$2
		given_value=$3
		project_id=$4
		numbers='^[0-9]+$'
		if [[ $given_value =~ $numbers ]] || [[ $given_value == "default" ]] ; then true; else 
			echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - found non-numerical value in "${master_file}" where expected to be numerical"
			echo -e "   Project: "${project_id}
			echo -e "   Column: "${master_column}
			echo -e "   Given value: "${given_value}
			echo -e ""${COL2}
			echo "error" >> $ERROR_MASTERPROJECT
		fi
	}
	check_value_num_in_master ${MASTERPROJECTS} "BIN" "${BIN}" "${PROJECT_ID}"
	check_value_num_in_master ${MASTERPROJECTS} "PADDING" "${PADDING}" "${PROJECT_ID}"

	#master project / segmentation
	# custom segmentation is defined in segmentation master
	if [[ ${SEGMENTATION_ID} != "default" ]]; then
		if [[ ${SEGMENTATION_ID} != "smart-"* ]]; then
			SEGMENTATION_HEADER=`cat ${SEGMENTATION} | head -1`
			if [[ ${SEGMENTATION_HEADER} =~ (^|[[:space:]])"${SEGMENTATION_ID}"($|[[:space:]]) ]]; then true; else
				echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - SEGMENTATION CONDITIONS: not recognized"
				echo -e "   Project: "${PROJECT_ID}
				echo -e "   Given SEGMENTATION_ID in the projects master ("${MASTERPROJECTS}"): "${SEGMENTATION_ID}
				echo -e "   This indicates usage of custom segmentation and requires values to be defined in the segmentation master ("${SEGMENTATION}") under a column called "${SEGMENTATION_ID}"."
				echo -e "   Segmentation conditions cannot be extracted as the column was not found. Please see README.txt for details on how to define segmentation."
				echo -e ""${COL2}
				echo "error" >> $ERROR_MASTERPROJECT
			fi
		fi
	fi




done

N_OF_ERRORS_MASTER=`cat ${ERROR_MASTERPROJECT} | wc -l`
rm ${ERROR_MASTERPROJECT}

if [[ ${N_OF_ERRORS_MASTER} == 0 ]]; then true; else
	echo -e ${COL1}${BEG}${COL2}
	echo -e ${COL1}"◼︎ $(date)"${COL2}
	echo -e ${COL1_WARN}"◼︎ CHECKPOINT ERROR - total number of "${N_OF_ERRORS_MASTER}" errors detected in the projects master "${MASTERPROJECTS}" (listed above)"
	echo -e "◼︎ Please resolve before processing."
	echo -e ${COL1}${END}${COL2}
	touch ./error
	exit
fi
#============================================================================




#CHECKPOINT - controls and samples master
#============================================================================

ERROR_MASTERCONTROLS="./error_controls"
touch $ERROR_MASTERCONTROLS
ERROR_MASTERSAMPLES="./error_samples"
touch $ERROR_MASTERSAMPLES

cat ${MASTERPROJECTS} | awk -F"\t" '$1=="yes"' | while read -r line || [[ -n "$line" ]]; do
	
	PROJECT_ID=`extract_from_master 'PROJECT_ID' ${MASTERPROJECTS}`
	PROJECT_TYPE=`extract_from_master 'PROJECT_TYPE' ${MASTERPROJECTS}`
	SEQ_TYPE=`extract_from_master 'SEQ_TYPE' ${MASTERPROJECTS}`
	CAPTURE_ID=`extract_from_master 'CAPTURE_ID' ${MASTERPROJECTS}`
	CAPTURE_FILE=`extract_from_master 'CAPTURE_FILE' ${MASTERPROJECTS}`
	GENOME_VERSION=`extract_from_master 'GENOME_VERSION' ${MASTERPROJECTS}`
	REF=`extract_from_master 'REF' ${MASTERPROJECTS}`
	SETTING_MODE=`extract_from_master 'SETTING_MODE' ${MASTERPROJECTS}`
	GENE_ANNOTATION=`extract_from_master 'GENE_ANNOTATION' ${MASTERPROJECTS}`
	BIN=`extract_from_master 'BIN' ${MASTERPROJECTS}`
	PADDING=`extract_from_master 'PADDING' ${MASTERPROJECTS}`
	WAY_TO_BAM=`extract_from_master 'WAY_TO_BAM' ${MASTERPROJECTS}`
	CTRL=`extract_from_master 'CTRL' ${MASTERPROJECTS}`
	CTRL_BAM_DIR=`extract_from_master 'CTRL_BAM_DIR' ${MASTERPROJECTS}`
	CTRL_BAM_PATTERN=`extract_from_master 'CTRL_BAM_PATTERN' ${MASTERPROJECTS}`
	CTRL_SEX_TEST=`extract_from_master 'CTRL_SEX_TEST' ${MASTERPROJECTS}`
	CTRL_PON_SEX_SELECT=`extract_from_master 'CTRL_PON_SEX_SELECT' ${MASTERPROJECTS}`
	SMPL_BAM_DIR=`extract_from_master 'SMPL_BAM_DIR' ${MASTERPROJECTS}`
	SMPL_BAM_PATTERN=`extract_from_master 'SMPL_BAM_PATTERN' ${MASTERPROJECTS}`
	SMPL_PON_SEX_SELECT=`extract_from_master 'SMPL_PON_SEX_SELECT' ${MASTERPROJECTS}`
	CN_FREQ_IN_CTRLS=`extract_from_master 'CN_FREQ_IN_CTRLS' ${MASTERPROJECTS}`
	GNOMAD_SELECTION=`extract_from_master 'GNOMAD_SELECTION' ${MASTERPROJECTS}`
	SEGMENTATION_ID=`extract_from_master 'SEGMENTATION_ID' ${MASTERPROJECTS}`
	SEGMENTATION_ID_USE=`extract_from_master 'SEGMENTATION_ID_USE' ${MASTERPROJECTS}`


	#function to check if provided value does not contain space character
	function check_space_in_master {
		master_file=$1
		master_column=$2
		given_value=$3
		project_id=$4
		main_id=$5
		case_id=$6
		error=$7
		if [[ "$given_value" = *" "* ]]; then
			echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - space character detected in "${master_file}
			echo -e "   Project: "${project_id}
			echo -e "   Column: "${master_column}
			echo -e "   Main id: "${main_id}
			echo -e "   Case id: "${case_id}
			echo -e "   Given value (in brackets to see all space characters): ("${given_value}")"
			echo -e ""${COL2}
			echo "error" >> $error
		fi
	}

	#function to check if provided value is not empty
	function check_empty_in_master_ctrl_smpl {
		master_file=$1
		master_column=$2
		given_value=$3
		project_id=$4
		main_id=$5
		case_id=$6
		error=$7
		if [[ -z "$given_value" ]]; then
			echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - unexpected empty value in "${master_file}
			echo -e "   Project: "${project_id}
			echo -e "   Column: "${master_column}
			echo -e "   Main id: "${main_id}
			echo -e "   Case id: "${case_id}
			echo -e ""${COL2}
			echo "error" >> $error
		fi
	}
	
	#function to check if provided value match expected values
	function check_value_in_master_ctrl_smpl {
		master_file=$1
		master_column=$2
		given_value=$3
		possible_values=$4
		project_id=$5
		main_id=$6
		case_id=$7
		error=$8
		if [[ ${possible_values} =~ (^|[[:space:]])"${given_value}"($|[[:space:]]) ]]; then true; else
				echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - unexpected value in "${master_file}
				echo -e "   Project: "${project_id}
				echo -e "   Column: "${master_column}
				echo -e "   Main id: "${main_id}
				echo -e "   Case id: "${case_id}
				echo -e "   Given value: "${given_value}
				echo -e "   Possible values: "${possible_values}" (see README.txt for details)"
				echo -e ""${COL2}
				echo "error" >> $error
		fi
	}



	#1) check master of controls
	#============================================================================
	#continue only when controls are expected
	if [[ ${CTRL} == "yes" ]]; then

		Nx=`cat ${MASTERCONTROLS} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION)print}' | wc -l`
		N=$(echo $Nx | tr -d ' ')

		#check if controls are filled when expected
		if [[ ${N} == 0 ]]; then
			echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - no controls found in the controls master "${MASTERCONTROLS}", but expected as the CTRL column in the projects master is yes"
			echo -e "   Project: "${PROJECT_ID}
			echo -e "   Please check if wanted controls have the first column (INCLUDE) as: yes"
			echo -e "   Please check if the following columns in the controls master are matching values defined in the projects master: PROJECT_ID, CAPTURE_ID, GENOME_VERSION"
			echo -e ""${COL2}
			echo "error" >> $ERROR_MASTERCONTROLS
		fi

		#continue to check the master file only when controls expected and filled
		if [[ ${N} != 0 ]]; then

		 	cat ${MASTERCONTROLS} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION)print}' | while read -r line || [[ -n "$line" ]]; do
				MAIN_ID=`extract_from_master 'MAIN_ID' ${MASTERCONTROLS}`
				CTRL_ID=`extract_from_master 'CTRL_ID' ${MASTERCONTROLS}`
				CTRL_SEX=`extract_from_master 'CTRL_SEX' ${MASTERCONTROLS}`
				CTRL_BAM_PATH=`extract_from_master 'CTRL_PATH_TO_BAM' ${MASTERCONTROLS}`

				#check space characters
				check_space_in_master "${MASTERCONTROLS}" "MAIN_ID" "${MAIN_ID}" "${PROJECT_ID}" "${MAIN_ID}" "${CTRL_ID}" "${ERROR_MASTERCONTROLS}"
				check_space_in_master "${MASTERCONTROLS}" "CTRL_ID" "${CTRL_ID}" "${PROJECT_ID}" "${MAIN_ID}" "${CTRL_ID}" "${ERROR_MASTERCONTROLS}"
				check_space_in_master "${MASTERCONTROLS}" "CTRL_SEX" "${CTRL_SEX}" "${PROJECT_ID}" "${MAIN_ID}" "${CTRL_ID}" "${ERROR_MASTERCONTROLS}"
				check_space_in_master "${MASTERCONTROLS}" "CTRL_BAM_PATH" "${CTRL_BAM_PATH}" "${PROJECT_ID}" "${MAIN_ID}" "${CTRL_ID}" "${ERROR_MASTERCONTROLS}"

				#check empty values
				check_empty_in_master_ctrl_smpl "${MASTERCONTROLS}" "MAIN_ID" "${MAIN_ID}" "${PROJECT_ID}" "${MAIN_ID}" "${CTRL_ID}" "${ERROR_MASTERCONTROLS}"
				check_empty_in_master_ctrl_smpl "${MASTERCONTROLS}" "CTRL_ID" "${CTRL_ID}" "${PROJECT_ID}" "${MAIN_ID}" "${CTRL_ID}" "${ERROR_MASTERCONTROLS}"

				#check expected values
				CTRL_SEX_VALUES=("M" "F" "unk")
				check_value_in_master_ctrl_smpl "${MASTERCONTROLS}" "CTRL_SEX" "${CTRL_SEX}" "${CTRL_SEX_VALUES[*]}" "${PROJECT_ID}" "${MAIN_ID}" "${CTRL_ID}" "${ERROR_MASTERCONTROLS}"

				#check bam files
				[[ ${WAY_TO_BAM} == "find_in_dir" ]] && BAMFILE=`find ${CTRL_BAM_DIR} -type f -name ${CTRL_ID}*${CTRL_BAM_PATTERN}`
				
				if [[ ${WAY_TO_BAM} == "find_in_dir" ]]; then
					BAMFILE_N=$(echo ${BAMFILE} | wc -w | tr -d ' ')
					if [[ $BAMFILE_N == 0 ]]; then
						echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - could not find BAM file filled in "${MASTERCONTROLS}
						echo -e "   Expected 1 BAM file but found "${BAMFILE_N}
						echo -e "   Project: "${PROJECT_ID}
						echo -e "   Main id: "${MAIN_ID}
						echo -e "   Case id: "${CTRL_ID}
						echo -e "   BAM file was indicated to be found in folder (defined in the projects master): "${CTRL_BAM_DIR}
						echo -e "   By pattern (defined in the controls master): "${CTRL_ID}" ... "${CTRL_BAM_PATTERN}
						echo -e ""${COL2}
						echo "error" >> $ERROR_MASTERCONTROLS
					fi
					if [[ $BAMFILE_N != 0 ]] && [[ $BAMFILE_N != 1 ]]; then
						echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - found more BAM files matching given pattern"
						echo -e "   Expected 1 BAM file but found "${BAMFILE_N}
						echo -e "   Project: "${PROJECT_ID}
						echo -e "   Main id: "${MAIN_ID}
						echo -e "   Case id: "${CTRL_ID}
						echo -e "   BAM file was indicated to be found in folder (defined in the projects master): "${CTRL_BAM_DIR}
						echo -e "   By pattern (defined in the controls master): "${CTRL_ID}" ... "${CTRL_BAM_PATTERN}
						echo -e ""${COL2}
						echo "error" >> $ERROR_MASTERCONTROLS
					fi
				fi
					
				[[ ${WAY_TO_BAM} == "absolute" ]] || [[ ${WAY_TO_BAM} == "default" ]] && BAMFILE=${CTRL_BAM_DIR}${CTRL_BAM_PATH}

				if [[ ${WAY_TO_BAM} == "absolute" ]] || [[ ${WAY_TO_BAM} == "default" ]]; then
					if ! ls ${BAMFILE} 1> /dev/null 2>&1; then
						echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - could not find BAM file filled in "${MASTERCONTROLS}
						echo -e "   Expected 1 BAM file but found 0"
						echo -e "   Project: "${PROJECT_ID}
						echo -e "   Main id: "${MAIN_ID}
						echo -e "   Case id: "${CTRL_ID}
						echo -e "   BAM file was indicated to be found in folder (defined in the projects master): "${CTRL_BAM_DIR}
						echo -e "   File name (defined in the controls master): "${CTRL_BAM_PATH}
						echo -e ""${COL2}
						echo "error" >> $ERROR_MASTERCONTROLS
					fi
				fi
			
				# echo $ERROR_MASTERCONTROLS > ./temp_ctrl
			done
		fi
	fi
	#1) check master of controls
	#============================================================================




	#2) check master of samples
	#============================================================================

	#check if at least one sample is filled
	Nx=`cat ${MASTERSAMPLES} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION)print}' | wc -l`
	N=$(echo $Nx | tr -d ' ')


	#check if samples are filled when expected
	if [[ ${N} == 0 ]]; then
		echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - no sample found in the samples master "${MASTERSAMPLES}" to be run"
		echo -e "   Project: "${PROJECT_ID}
		echo -e "   Please check if wanted samples have the first column (INCLUDE) as: yes"
		echo -e "   Please check if the following columns in the samples master are matching values defined in the projects master: PROJECT_ID, CAPTURE_ID, GENOME_VERSION"
		echo -e ""${COL2}
		echo "error" >> $ERROR_MASTERSAMPLES
	fi


	#continue to check the master file only when at least one sample is recognized
	if [[ ${N} != 0 ]]; then

	 	cat ${MASTERSAMPLES} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION)print}' | while read -r line || [[ -n "$line" ]]; do
			
		 	MAIN_ID=`extract_from_master 'MAIN_ID' ${MASTERSAMPLES}`
		 	SAMPLE_1_ID=`extract_from_master 'SAMPLE1_ID' ${MASTERSAMPLES}`
		 	SAMPLE_1_TYPE=`extract_from_master 'SAMPLE1_TYPE' ${MASTERSAMPLES}`
		 	SAMPLE_1_SEX=`extract_from_master 'SAMPLE1_SEX' ${MASTERSAMPLES}`
		 	SAMPLE_1_BAM_PATH=`extract_from_master 'SAMPLE1_PATH_TO_BAM' ${MASTERSAMPLES}`

			#check space characters
			check_space_in_master "${MASTERSAMPLES}" "MAIN_ID" "${MAIN_ID}" "${PROJECT_ID}" "${MAIN_ID}" "${SAMPLE_1_ID}" "${ERROR_MASTERSAMPLES}"
			check_space_in_master "${MASTERSAMPLES}" "SAMPLE_1_ID" "${SAMPLE_1_ID}" "${PROJECT_ID}" "${MAIN_ID}" "${SAMPLE_1_ID}" "${ERROR_MASTERSAMPLES}"
			check_space_in_master "${MASTERSAMPLES}" "SAMPLE_1_TYPE" "${SAMPLE_1_TYPE}" "${PROJECT_ID}" "${MAIN_ID}" "${SAMPLE_1_ID}" "${ERROR_MASTERSAMPLES}"
			check_space_in_master "${MASTERSAMPLES}" "SAMPLE_1_SEX" "${SAMPLE_1_SEX}" "${PROJECT_ID}" "${MAIN_ID}" "${SAMPLE_1_ID}" "${ERROR_MASTERSAMPLES}"
			check_space_in_master "${MASTERSAMPLES}" "SAMPLE_1_BAM_PATH" "${SAMPLE_1_BAM_PATH}" "${PROJECT_ID}" "${MAIN_ID}" "${SAMPLE_1_ID}" "${ERROR_MASTERSAMPLES}"

			#check empty values
			check_empty_in_master_ctrl_smpl ${MASTERSAMPLES} "MAIN_ID" "${MAIN_ID}" "${PROJECT_ID}" "${MAIN_ID}" "${SAMPLE_1_ID}" "${ERROR_MASTERSAMPLES}"
			check_empty_in_master_ctrl_smpl ${MASTERSAMPLES} "SAMPLE_1_ID" "${SAMPLE_1_ID}" "${PROJECT_ID}" "${MAIN_ID}" "${SAMPLE_1_ID}" "${ERROR_MASTERSAMPLES}"
			check_empty_in_master_ctrl_smpl ${MASTERSAMPLES} "SAMPLE_1_TYPE" "${SAMPLE_1_TYPE}" "${PROJECT_ID}" "${MAIN_ID}" "${SAMPLE_1_ID}" "${ERROR_MASTERSAMPLES}"

			#check expected values
			SMPL_SEX_VALUES=("M" "F" "unk")
			check_value_in_master_ctrl_smpl ${MASTERSAMPLES} "SAMPLE_1_SEX" "${SAMPLE_1_SEX}" "${SMPL_SEX_VALUES[*]}" "${PROJECT_ID}" "${MAIN_ID}" "${SAMPLE_1_ID}" "${ERROR_MASTERSAMPLES}"


			#check bam files
			[[ ${WAY_TO_BAM} == "find_in_dir" ]] && BAMFILE=`find ${SMPL_BAM_DIR} -type f -name ${SAMPLE_1_ID}*${SMPL_BAM_PATTERN}`
			
			if [[ ${WAY_TO_BAM} == "find_in_dir" ]]; then
				BAMFILE_N=$(echo ${BAMFILE} | wc -w | tr -d ' ')
				if [[ $BAMFILE_N == 0 ]]; then
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - could not find BAM file filled in "${MASTERSAMPLES}
					echo -e "   Expected 1 BAM file but found "${BAMFILE_N}
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Main id: "${MAIN_ID}
					echo -e "   Case id: "${SAMPLE_1_ID}
					echo -e "   BAM file was indicated to be found in folder (defined in the projects master): "${SMPL_BAM_DIR}
					echo -e "   By pattern (defined in the samples master): "${SAMPLE_1_ID}" ... "${SMPL_BAM_PATTERN}
					echo -e ""${COL2}
					echo "error" >> $ERROR_MASTERSAMPLES
				fi
				if [[ $BAMFILE_N != 0 ]] && [[ $BAMFILE_N != 1 ]]; then
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - found more BAM files matching given pattern"${MASTERSAMPLES}
					echo -e "   Expected 1 file but found "${BAMFILE_N}
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Main id: "${MAIN_ID}
					echo -e "   Case id: "${SAMPLE_1_ID}
					echo -e "   BAM file was indicated to be found in folder (defined in the projects master): "${CTRL_BAM_DIR}
					echo -e "   By pattern (defined in the samples master): "${SAMPLE_1_ID}" ... "${SMPL_BAM_PATTERN}
					echo -e ""${COL2}
					echo "error" >> $ERROR_MASTERSAMPLES
				fi
			fi
				
			[[ ${WAY_TO_BAM} == "absolute" ]] || [[ ${WAY_TO_BAM} == "default" ]] && BAMFILE=${SMPL_BAM_DIR}${SAMPLE_1_BAM_PATH}

			if [[ ${WAY_TO_BAM} == "absolute" ]] || [[ ${WAY_TO_BAM} == "default" ]]; then
				if ! ls ${BAMFILE} 1> /dev/null 2>&1; then
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - could not find BAM file filled in "${MASTERSAMPLES}
					echo -e "   Expected 1 BAM file but found 0"
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Main id: "${MAIN_ID}
					echo -e "   Case id: "${SAMPLE_1_ID}
					echo -e "   BAM file was indicated to be found in folder (defined in the projects master): "${SMPL_BAM_DIR}
					echo -e "   File name (defined in the samples master): "${SAMPLE_1_BAM_PATH}
					echo -e ""${COL2}
					echo "error" >> $ERROR_MASTERSAMPLES
				fi
			fi		
		done
	fi
	#2) check master of samples
	#============================================================================

	#3) check capture file
	#============================================================================
	# part of run_01_prep.sh
	#============================================================================

	#4) check sex master
	#============================================================================
	# part of run_01_prep.sh
	#============================================================================

	#5) check segmentation master - numerical values expected
	#============================================================================
	# part of run_01_prep.sh 
	#============================================================================

done

N_OF_ERRORS_CONTROLS=`cat ${ERROR_MASTERCONTROLS} | wc -l`
rm ${ERROR_MASTERCONTROLS}
N_OF_ERRORS_SAMPLES=`cat ${ERROR_MASTERSAMPLES} | wc -l`
rm ${ERROR_MASTERSAMPLES}

#change code below for ctrl and smpl errors
if [[ ${N_OF_ERRORS_CONTROLS} != 0 ]] || [[ ${N_OF_ERRORS_SAMPLES} != 0 ]]; then
	echo -e ${COL1}${BEG}${COL2}
	echo -e ${COL1}"◼︎ $(date)"${COL2}
	echo -e ${COL1_WARN}"◼︎ CHECKPOINT ERROR - errors detected in the controls master and/or the samples master (listed above)"
	echo -e ${COL1_WARN}"◼︎ Number of errors in the controls master "${MASTERCONTROLS}": "${N_OF_ERRORS_CONTROLS}
	echo -e ${COL1_WARN}"◼︎ Number of errors in the samples master "${MASTERSAMPLES}": "${N_OF_ERRORS_SAMPLES}	
	echo -e "◼︎ Please resolve before processing."
	echo -e ${COL1}${END}${COL2}
	touch ./error
	exit
fi

echo -e ${COL1}${BEG}${COL2}
echo -e ${COL1}"◼︎ $(date)"${COL2}
echo -e ${COL1}"◼︎ Initial testing completed with no error message"${COL2}
echo -e ${COL1}${END}${COL2}



#============================================================================