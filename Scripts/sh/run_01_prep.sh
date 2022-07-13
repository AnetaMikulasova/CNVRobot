#!/bin/bash

source ${1}

echo -e ${COL1}${BEG}${COL2}
echo -e ${COL1}"◼︎ $(date)"${COL2}
echo -e ${COL1}"◼︎ Controls preparation ..."${COL2}
echo -e ${COL1}${END}${COL2}

##### preparation of intervals, controls and gnomad database
#####=======================================================================================================

#LOOP - LINE oF PROJECTS
#============================================================================
cat ${MASTERPROJECTS} | awk -F"\t" '$1=="yes"' | while read -r line || [[ -n "$line" ]]; do
	
	#prevent looping if error occured during previous round
	if [[ ${CHECKPOINTS} == "yes" ]]; then 
		if ls ./error 1> /dev/null 2>&1; then exit; fi
	fi

	PROJECT_ID=`extract_from_master 'PROJECT_ID' ${MASTERPROJECTS}`
	PROJECT_TYPE=`extract_from_master 'PROJECT_TYPE' ${MASTERPROJECTS}`
	SEQ_TYPE=`extract_from_master 'SEQ_TYPE' ${MASTERPROJECTS}`
	CAPTURE_ID=`extract_from_master 'CAPTURE_ID' ${MASTERPROJECTS}`
	CAPTURE_FILE=`extract_from_master 'CAPTURE_FILE' ${MASTERPROJECTS}`
	GENOME_VERSION=`extract_from_master 'GENOME_VERSION' ${MASTERPROJECTS}`
	REF=`extract_from_master 'REF' ${MASTERPROJECTS}`
	READ_COUNT_MODE=`extract_from_master 'READ_COUNT_MODE' ${MASTERPROJECTS}`
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

	echo -e ${COL1}${BEG}${COL2}
	echo -e ${COL1}"◼︎ $(date)"${COL2}
	echo -e ${COL1}"◼︎ PREPROCESSING - project: "${PROJECT_ID}" - segmentation: "${SEGMENTATION_ID}${COL2}
	echo -e ${COL1}${END}${COL2}

	###define default setting
	###====================================================================================
	[ ${READ_COUNT_MODE} == "default" ] && READ_COUNT_MODE="F"
	[ ${WAY_TO_BAM} == "default" ] && WAY_TO_BAM="absolute"
	[ ${CTRL_BAM_PATTERN} == "default" ] && CTRL_BAM_PATTERN="na"
	[ ${CTRL_PON_SEX_SELECT} == "default" ] && CTRL_PON_SEX_SELECT="matched"
	[ ${SMPL_BAM_PATTERN} == "default" ] && SMPL_BAM_PATTERN="na"
	[ ${SMPL_PON_SEX_SELECT} == "default" ] && SMPL_PON_SEX_SELECT="matched_each"
	[ ${SEGMENTATION_ID_USE} == "default" ] && SEGMENTATION_ID_USE="segm_id"

	#define segmentation condition based on PoN selected for samples
	[ ${SEGMENTATION_ID} == "default" ] && [ ${SMPL_PON_SEX_SELECT} == "mixed" ] && SEGMENTATION_ID="smart-0.85"
	[ ${SEGMENTATION_ID} == "default" ] && [ ${SMPL_PON_SEX_SELECT} != "mixed" ] && SEGMENTATION_ID="smart-0.65"

	#define gnomad selection based on sequencing type
	[ ${SEQ_TYPE} == "WGS" ] && [ ${GNOMAD_SELECTION} == "default" ] && GNOMAD_SELECTION="af_filter_and_ctrl_denois"
	[ ${SEQ_TYPE} != "WGS" ] && [ ${GNOMAD_SELECTION} == "default" ] && GNOMAD_SELECTION="capture_filter_and_ctrl_denois"

	#define bin and pad based on sequencing type
	#WGS
	[ ${SEQ_TYPE} == "WGS" ] && [ ${BIN} == "default" ] && BIN=1000
	[ ${SEQ_TYPE} == "WGS" ] && [ ${PADDING} == "default" ] && PADDING=0
	#WES or TS
	if [[ ${SEQ_TYPE} != "WGS" ]]; then
		if [[ ${BIN} == "default" ]]; then
			BIN=1000
		fi
		if [[ ${PADDING} == "default" ]]; then
			#to calculate padding, enter first control and calculate read length
			line=`cat ${MASTERCONTROLS} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION)print}' | head -1`
			MAIN_ID=`extract_from_master 'MAIN_ID' ${MASTERCONTROLS}`
			CTRL_ID=`extract_from_master 'CTRL_ID' ${MASTERCONTROLS}`
			CTRL_SEX=`extract_from_master 'CTRL_SEX' ${MASTERCONTROLS}`
			CTRL_BAM_PATH=`extract_from_master 'CTRL_PATH_TO_BAM' ${MASTERCONTROLS}`
			[ ${WAY_TO_BAM} == "find_in_dir" ] && BAMFILE=`find ${CTRL_BAM_DIR} -type f -name ${CTRL_ID}*${CTRL_BAM_PATTERN}`
			[ ${WAY_TO_BAM} == "absolute" ] && BAMFILE=${CTRL_BAM_DIR}${CTRL_BAM_PATH}
			PADDING=`samtools view $BAMFILE | awk '{print length($10)}' | head -10000 | sort -nr | head -1`
		fi
	fi
	###====================================================================================
	
	###extract mode
	###====================================================================================
	READ_COUNT_FILTER=$(echo $READ_COUNT_MODE | awk -F"-" '{print $1}')
	GERMLINE_MODE=$(echo $READ_COUNT_MODE | awk -F"-" '{print $2}')
	QUICK_MODE=$(echo $READ_COUNT_MODE | awk -F"-" '{print $3}')
	QUICK_MODE_PLOT=$(echo $READ_COUNT_MODE | awk -F"-" '{print $4}')
	if [[ ${GERMLINE_MODE} != "G" ]]; then 
		MODE_PATTERN=$(echo $READ_COUNT_FILTER)
	fi
	if [[ ${GERMLINE_MODE} == "G" ]]; then 
		MODE_PATTERN=$(echo $READ_COUNT_FILTER"-"$GERMLINE_MODE)
	fi
	###====================================================================================

	#count number of controls
	Nx=`cat ${MASTERCONTROLS} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION)print}' | wc -l`
	[ ${CTRL} == "yes" ] && N=$(echo $Nx | tr -d ' ')
	[ ${CTRL} == "no" ] && N=0

	SUPPORT_FILES=${SUPPORT_FILES_DIR}"/"${PROJECT_ID}_${CAPTURE_ID}"_"${GENOME_VERSION}"_"${MODE_PATTERN}"_bin"${BIN}"bp_pad"${PADDING}"bp/"

	#assigned folders	
	INTERVALS=${SUPPORT_FILES}"Capture_Intervals/"
	CTRL_SEX_DIR=${SUPPORT_FILES}"Controls/sex/"
	CTRL_COUNT_HDF_DIR=${SUPPORT_FILES}"Controls/counts/hdf5/"
	CTRL_COUNT_TSV_DIR=${SUPPORT_FILES}"Controls/counts/tsv/"
	CTRL_MAF_DIR=${SUPPORT_FILES}"Controls/maf/"
	CTRL_PoN_DIR=${SUPPORT_FILES}"Controls/PoN/NofCTRL="${N}"/"
	CTRL_STANDARD_DIR=${SUPPORT_FILES}"Controls/standardized/NofCTRL="${N}"/"
	CTRL_DENOIS_DIR=${SUPPORT_FILES}"Controls/denois/NofCTRL="${N}"/"
	CTRL_SEGMENT_DIR=${SUPPORT_FILES}"Controls/segmentation/NofCTRL="${N}"/"
	CTRL_QC_DIR=${SUPPORT_FILES}"Controls/qc/NofCTRL="${N}"/"
	CTRL_CN_DIR=${SUPPORT_FILES}"Controls/CN/"

	CONTIGS=${INTERVALS}${GENOME_VERSION}"_contigs.bed"
	CONTIGS_PROCESSED=${INTERVALS}${GENOME_VERSION}"_contigs.interval_list"
	CONTIGS_SIZES=${INTERVALS}${GENOME_VERSION}"_contigs_sizes.txt"
	# CAPTURE_FILE_PROCESSED=${INTERVALS}${CAPTURE_ID}"_sorted.interval_list"
	CAPTURE_FILE_FILTERED=${INTERVALS}${CAPTURE_ID}"_filtered.bed"
	CAPTURE_FILE_SORTED=${INTERVALS}${CAPTURE_ID}"_sorted.interval_list"

	# INTERVALS_PROCESSED_TEMP=${INTERVALS}${CAPTURE_ID}"_preprocessed_intervals_temp.interval_list"
	INTERVALS_PROCESSED=${INTERVALS}${CAPTURE_ID}"_preprocessed_intervals.interval_list"
	INTERVALS_ANNOTATED=${INTERVALS}${CAPTURE_ID}"_annotated_intervals.interval_list"

	#reference indexes
	REF_BASE="${REF%.*}"
	REF_DICT=${REF_BASE}".dict"
	REF_FAIDX=${REF}".fai"
	

	#recognize if user gives segmentation
	if [[ ${SEGMENTATION_ID} != "smart-"* ]]; then
		
		# segmentation - process conditions
		echo -e ${COL1}${BEG}${COL2}
		echo -e ${COL1}"◼︎ $(date)"${COL2}
		echo -e ${COL1}"◼︎ PREPROCESSING - preparing segmentation conditions"${COL2}
		echo -e ${COL1}${END}${COL2}
		SEGMENTATION_TEMP=${SEGMENTATION}"_temp.txt"
		if ls ${SEGMENTATION_TEMP} 1> /dev/null 2>&1; then
			rm ${SEGMENTATION_TEMP}
		fi
		Rscript --vanilla ${R}"segmentation_conditions.R" \
				${SEGMENTATION} \
				${SEGMENTATION_ID} \
				${SEGMENTATION_TEMP}
		# SEGMENTATION=${SEGMENTATION_TEMP}

		# segmentation - get values
		NOGAPS=250000000
		SEGM_DIFFERENCE=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_DIFFERENCE")print $2}'`
		SEGM_MINSIZE=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_MINSIZE")print $2}'`
		SEGM_MINSIZE_SUB=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_MINSIZE_SUB")print $2}'`
		[ ${SEGM_MINSIZE_SUB} == "none" ] && SEGM_MINSIZE_SUB=${SEGM_MINSIZE}
		SEGM_MINPROBE=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_MINPROBE")print $2}'`
		SEGM_MINPROBE_SUB=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_MINPROBE_SUB")print $2}'`
		[ ${SEGM_MINPROBE_SUB} == "none" ] && SEGM_MINPROBE_SUB=${SEGM_MINPROBE}
		SEGM_MINKEEP=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_MINKEEP")print $2}'`
		SEGM_MINLOSS=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_MINLOSS")print $2}'`
		SEGM_MINLOSS_SUB=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_MINLOSS_SUB")print $2}'`
		[ ${SEGM_MINLOSS_SUB} == "none" ] && SEGM_MINLOSS_SUB=${SEGM_MINLOSS}
		SEGM_MINLOSS_BIAL=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_MINLOSS_BIAL")print $2}'`
		SEGM_MINGAIN=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_MINGAIN")print $2}'`
		SEGM_MINGAIN_SUB=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_MINGAIN_SUB")print $2}'`
		[ ${SEGM_MINGAIN_SUB} == "none" ] && SEGM_MINGAIN_SUB=${SEGM_MINGAIN}
		SEGM_GAP=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_GAP")print $2}'`
		[ ${SEGM_GAP} == "none" ] && SEGM_GAP=${NOGAPS}
		SEGM_SMOOTHPERC=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_SMOOTHPERC")print $2}'`
		SEGM_AFDIF=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_AFDIF")print $2}'`
		SEGM_AFSIZE=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_AFSIZE")print $2}'`
		SEGM_AFPROBE=`cat ${SEGMENTATION_TEMP} | awk -F"\t" '{if ($1=="SEGM_AFPROBE")print $2}'`
		

		SEGMENTATION_FULL=$SEGM_DIFFERENCE$x$SEGM_MINSIZE$x$SEGM_MINSIZE_SUB$x$SEGM_MINPROBE$x$SEGM_MINPROBE_SUB$x$SEGM_MINKEEP$x$SEGM_MINLOSS$x$SEGM_MINLOSS_SUB$x$SEGM_MINLOSS_BIAL$x$SEGM_MINGAIN$x$SEGM_MINGAIN_SUB$x$SEGM_GAP$x$SEGM_SMOOTHPERC$x$SEGM_AFDIF$x$SEGM_AFSIZE$x$SEGM_AFPROBE
		
		[ ${SEGMENTATION_ID_USE} == "full" ] && SEGMENTATIONCONDITIONS=${SEGMENTATION_FULL}
		[ ${SEGMENTATION_ID_USE} == "segm_id" ] && SEGMENTATIONCONDITIONS=${SEGMENTATION_ID}
	

		#CHECKPOINT - segmentation conditions 
		#============================================================================
		if [[ ${CHECKPOINTS} == "yes" ]]; then
			#function to check if values values are numerical
			function check_segm_value_num {
				segmentation_condition=$1
				given_value=$2
				numbers='^[+-]?[0-9]+([.][0-9]+)?$'
				if [[ $given_value =~ $numbers ]]; then true; else 
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - SEGMENTATION CONDITIONS: some condition with unexpected value"
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Segmentation conditions were indicated to be defined as custom values in the segmentation master ("${SEGMENTATION}") under column called "${SEGMENTATION_ID}"."
					echo -e "   Some segmentation condition was not given correctly (missing or not numerical value)."
					echo -e "   Segmentation condition: "${segmentation_condition}
					echo -e "   Extracted value: "${given_value}
					echo -e "   Please check this condition in the segmentation master ("${SEGMENTATION}")."
					echo -e ""${COL2}
					touch ./error
					exit
				fi
			}
			check_segm_value_num "SEGM_DIFFERENCE" "${SEGM_DIFFERENCE}"
			check_segm_value_num "SEGM_MINSIZE" "${SEGM_MINSIZE}"
			check_segm_value_num "SEGM_MINSIZE_SUB" "${SEGM_MINSIZE_SUB}"
			check_segm_value_num "SEGM_MINPROBE" "${SEGM_MINPROBE}"
			check_segm_value_num "SEGM_MINPROBE_SUB" "${SEGM_MINPROBE_SUB}"
			check_segm_value_num "SEGM_MINKEEP" "${SEGM_MINKEEP}"
			check_segm_value_num "SEGM_MINLOSS" "${SEGM_MINLOSS}"
			check_segm_value_num "SEGM_MINLOSS_SUB" "${SEGM_MINLOSS_SUB}"
			check_segm_value_num "SEGM_MINLOSS_BIAL" "${SEGM_MINLOSS_BIAL}"
			check_segm_value_num "SEGM_MINGAIN" "${SEGM_MINGAIN}"
			check_segm_value_num "SEGM_MINGAIN_SUB" "${SEGM_MINGAIN_SUB}"
			check_segm_value_num "SEGM_GAP" "${SEGM_GAP}"
			check_segm_value_num "SEGM_SMOOTHPERC" "${SEGM_SMOOTHPERC}"
			check_segm_value_num "SEGM_AFDIF" "${SEGM_AFDIF}"
			check_segm_value_num "SEGM_AFSIZE" "${SEGM_AFSIZE}"
			check_segm_value_num "SEGM_AFPROBE" "${SEGM_AFPROBE}"
		fi
		#CHECKPOINT - segmentation conditions 
		#============================================================================

	#recognize if user gives segmentation
	fi




	# ### prepare reference contigs and bed files
	# ###====================================================================================

		### create dict file for the reference - needed for sorting bed files
		#============================================================================
		if ! ls ${REF_DICT} 1> /dev/null 2>&1; then
			echo -e ${COL1}${BEG}${COL2}
			echo -e ${COL1}"◼︎ $(date)"${COL2}
			echo -e ${COL1}"◼︎ PREPROCESSING - preparing dict file for fasta reference file"${COL2}
			echo -e ${COL1}${END}${COL2}
			${PICARD} CreateSequenceDictionary \
				R=${REF} \
				O=${REF_DICT}
 			if [[ ${CHECKPOINTS} == "yes" ]]; then 
				if ! ls ${REF_DICT} 1> /dev/null 2>&1; then
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - REFERENCE PREP: file not generated"
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Stage: indexing reference"
					echo -e "   File not generated: "${REF_DICT}
					echo -e ""${COL2}
					touch ./error
					exit
				fi
			fi
		 fi
		#============================================================================

		### index the reference - needed for sorting bed files
		#============================================================================
		if ! ls ${REF_FAIDX} 1> /dev/null 2>&1; then
			echo -e ${COL1}${BEG}${COL2}
			echo -e ${COL1}"◼︎ $(date)"${COL2}
			echo -e ${COL1}"◼︎ PREPROCESSING - preparing fasta index file"${COL2}
			echo -e ${COL1}${END}${COL2}
			${SAMTOOLS} faidx ${REF}
 			if [[ ${CHECKPOINTS} == "yes" ]]; then 
				if ! ls ${REF_FAIDX} 1> /dev/null 2>&1; then
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - REFERENCE PREP: file not generated"
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Stage: indexing reference"
					echo -e "   File not generated: "${REF_FAIDX}
					echo -e ""${COL2}
					touch ./error
					exit
				fi
			fi
		 fi
		# ============================================================================

		### export contigs from ref fai index
		#============================================================================
		if ! ls ${CONTIGS} 1> /dev/null 2>&1; then
			echo -e ${COL1}${BEG}${COL2}
			echo -e ${COL1}"◼︎ $(date)"${COL2}
			echo -e ${COL1}"◼︎ PREPROCESSING - exporting contigs and their sizes from reference"${COL2}
			echo -e ${COL1}${END}${COL2}
			mkdir -p ${INTERVALS}
			Rscript --vanilla ${R}contigs.R \
				${REF_FAIDX} \
				${CONTIGS} \
				${CONTIGS_SIZES}
			if [[ ${CHECKPOINTS} == "yes" ]]; then
				if ls "status_ok" 1> /dev/null 2>&1; then
					rm "status_ok"
				else
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: R script was halted"
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Stage: extraction of contigs from reference"
					echo -e "   File not generated: "${CONTIGS} " and/or "${CONTIGS_SIZES}
					echo -e "   R script: contigs.R"
					echo -e ""${COL2}
					touch ./error
					exit
				fi
				if ! ls ${CONTIGS} 1> /dev/null 2>&1; then
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Stage: extraction of contigs from reference"
					echo -e "   File not generated: "${CONTIGS}
					echo -e ""${COL2}
					touch ./error
					exit
				fi
				if ! ls ${CONTIGS_SIZES} 1> /dev/null 2>&1; then
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Stage: extraction of contigs from reference"
					echo -e "   File not generated: "${CONTIGS_SIZES}
					echo -e ""${COL2}
					touch ./error
					exit
				fi
			fi
		fi
		#============================================================================

	# ###====================================================================================


	### prepare intervals for specific capture by given bin and padding, and annotate for %GC
	###====================================================================================

		### process intervals by given bin and padding
		#============================================================================
		if ! ls ${INTERVALS_PROCESSED} 1> /dev/null 2>&1; then

			#a) WES or any other targeted sequencing with given capture file
			if [ ${SEQ_TYPE} != "WGS" ]; then
				#filter capture intervals for required contigs only
				echo -e ${COL1}${BEG}${COL2}
				echo -e ${COL1}"◼︎ $(date)"${COL2}
				echo -e ${COL1}"◼︎ PREPROCESSING - filtering capture file for required contigs only"${COL2}
				echo -e ${COL1}${END}${COL2}
				dos2unix ${CAPTURE_FILE}
				${BEDTOOLS} intersect -a ${CAPTURE_FILE} -b ${CONTIGS} > ${CAPTURE_FILE_FILTERED}

				#CHECKPOINT - capture file
				#============================================================================
				if [[ ${CHECKPOINTS} == "yes" ]]; then
					#check three columns in capture file
					NofCOLUMNS=`cat ${CAPTURE_FILE} | awk -F"\t" '{print NF; exit}'`
					if [[ "$NofCOLUMNS" -lt 3 ]]; then
						echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CAPTURE FILE: found less than 3 columns in capture BED file"
						echo -e "   Project: "${PROJECT_ID}
						echo -e "   Please check that capture file has 3 columns (contig, start and end) and no header."
						echo -e ""${COL2}
						touch ./error
						exit
					fi
					#checkpoint if at least some contigs are present in capture file after excluding those not present in reference (can identify situation when chr1, chr2... versus 1, 2...)
					NofLINES=`cat ${CAPTURE_FILE_FILTERED} | wc -l`
					CONTIG_in_REF=`cat ${REF_DICT} | grep "@SQ" | awk -F"\t|:" '{print $3}'`
					CONTIG_in_CAPTURE=`cat ${CAPTURE_FILE} | awk -F"\t" '{print $1}' | uniq`
					if [[ $NofLINES == 0 ]]; then
						echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CAPTURE FILE: no contigs in the capture file recognized within the reference file"
						echo -e "   Project: "${PROJECT_ID}
						echo -e "   Contigs in capture file: "${CONTIG_in_CAPTURE}
						echo -e "   Contigs in reference: "${CONTIG_in_REF}
						echo -e "   This often happens due to chr1, chr2, chr3 etc. versus 1, 2, 3 etc."
						echo -e "   Make sure that the capture file is in the same genome assembly as the reference (and all alignments)."
						echo -e ""${COL2}
						touch ./error
						exit
					fi
				fi
				#CHECKPOINT - capture file
				#============================================================================

				#sort capture file based on reference dict and transform from bed to interval.list
				echo -e ${COL1}${BEG}${COL2}
				echo -e ${COL1}"◼︎ $(date)"${COL2}
				echo -e ${COL1}"◼︎ PREPROCESSING - sorting capture file based on reference"${COL2}
				echo -e ${COL1}${END}${COL2}				
				${GATK} BedToIntervalList \
					-I ${CAPTURE_FILE_FILTERED} \
					-SD ${REF_DICT} \
					-O ${CAPTURE_FILE_SORTED}
				#process capture intervals by bin and padding
				${GATK} PreprocessIntervals \
					-R ${REF} \
					-L ${CAPTURE_FILE_SORTED} \
					--bin-length ${BIN} \
					--padding ${PADDING} \
					--interval-merging-rule OVERLAPPING_ONLY \
					-O ${INTERVALS_PROCESSED}
				rm ${CAPTURE_FILE_FILTERED}
				rm ${CAPTURE_FILE_SORTED}
			fi

			#b) WGS
			if [ ${SEQ_TYPE} == "WGS" ]; then
				# process contigs to interal.list
				echo -e ${COL1}${BEG}${COL2}
				echo -e ${COL1}"◼︎ $(date)"${COL2}
				echo -e ${COL1}"◼︎ PREPROCESSING - sorting assembled chromosomes file based on reference"${COL2}
				echo -e ${COL1}${END}${COL2}
				mkdir -p ${INTERVALS}
				${GATK} BedToIntervalList \
					-I ${CONTIGS} \
					-SD ${REF_DICT} \
					-O ${CONTIGS_PROCESSED}
				# process intervals by bin and padding
				echo -e ${COL1}${BEG}${COL2}
				echo -e ${COL1}"◼︎ $(date)"${COL2}
				echo -e ${COL1}"◼︎ "${CAPTURE_ID}" capture - processing intervals"${COL2}
				echo -e ${COL1}${END}${COL2}
				mkdir -p ${INTERVALS}
				${GATK} PreprocessIntervals \
					-R ${REF} \
					-L ${CONTIGS_PROCESSED} \
					--bin-length ${BIN} \
					--padding ${PADDING} \
					--interval-merging-rule OVERLAPPING_ONLY \
					-O ${INTERVALS_PROCESSED}
				rm ${CONTIGS_PROCESSED}
			fi
		
 			if [[ ${CHECKPOINTS} == "yes" ]]; then 
				if ! ls ${INTERVALS_PROCESSED} 1> /dev/null 2>&1; then
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - INTERVALS PREP: file not generated"
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Stage: preparing intervals"
					echo -e "   File not generated: "${INTERVALS_PROCESSED}
					echo -e ""${COL2}
					touch ./error
					exit
				fi
			fi

		fi
		#============================================================================

		### annotate intervals for GC correction
		#============================================================================
		if ! ls ${INTERVALS_ANNOTATED} 1> /dev/null 2>&1; then
			echo -e ${COL1}${BEG}${COL2}
			echo -e ${COL1}"◼︎ $(date)"${COL2}
			echo -e ${COL1}"◼︎ "${CAPTURE_ID}" capture - annotating intervals"${COL2}
			echo -e ${COL1}${END}${COL2}
			${GATK} AnnotateIntervals \
				-R ${REF} \
				-L ${INTERVALS_PROCESSED} \
				--interval-merging-rule OVERLAPPING_ONLY \
				-O ${INTERVALS_ANNOTATED}
 			if [[ ${CHECKPOINTS} == "yes" ]]; then 
				if ! ls ${INTERVALS_ANNOTATED} 1> /dev/null 2>&1; then
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - INTERVALS PREP: file not generated"
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Stage: annotation intervals"
					echo -e "   File not generated: "${INTERVALS_ANNOTATED}
					echo -e ""${COL2}
					touch ./error
					exit
				fi
			fi
		fi
		#============================================================================

	###====================================================================================

	### calculate segmentation conditions when smart option
	###====================================================================================
	# not earlier as it requires annotated intervals
	# SEGM_DIFFERENCE and SEGMENTATION_FULL generated later after data denoising
	if [[ ${SEGMENTATION_ID} == "smart-"* ]]; then
		SMART_SEGMENTATION=`Rscript --vanilla ${R}"segmentation_conditions_smart.R" \
			${SEGMENTATION_ID} \
			"na" \
			${INTERVALS_ANNOTATED} \
			${SEQ_TYPE}`

		#get rid off the R format
		SMART_SEGMENTATION=$(echo ${SMART_SEGMENTATION} | awk '{print $2}')
		SMART_SEGMENTATION=$(echo ${SMART_SEGMENTATION} | tr -d '\"')
		#extract values
		# SEGM_DIFFERENCE=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $1}') 
		SEGM_MINSIZE=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $2}') 
		SEGM_MINSIZE_SUB=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $3}') 
		SEGM_MINPROBE=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $4}') 
		SEGM_MINPROBE_SUB=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $5}') 
		SEGM_MINKEEP=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $6}') 
		SEGM_MINLOSS=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $7}') 
		SEGM_MINLOSS_SUB=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $8}') 
		SEGM_MINLOSS_BIAL=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $9}') 
		SEGM_MINGAIN=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $10}') 
		SEGM_MINGAIN_SUB=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $11}') 
		SEGM_GAP=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $12}') 
		SEGM_SMOOTHPERC=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $13}') 
		SEGM_AFDIF=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $14}') 
		SEGM_AFSIZE=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $15}') 
		SEGM_AFPROBE=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $16}') 

		#note: if segmenation is smart, segmentation contidions can't be "full" as it wouldn't be possible to find for report
		SEGMENTATIONCONDITIONS="SEGM-"${SEGMENTATION_ID}

		#CHECKPOINT - segmentation conditions 
		#============================================================================
		if [[ ${CHECKPOINTS} == "yes" ]]; then
			#function to check if values values are numerical
			function check_segm_value_num {
				segmentation_condition=$1
				given_value=$2
				numbers='^[+-]?[0-9]+([.][0-9]+)?$'
				if [[ $given_value =~ $numbers ]]; then true; else 
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - SEGMENTATION CONDITIONS: some condition with unexpected value"
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Segmentation conditions were indicated to be defined under smart segmentation, but format likely not provided correctly."
					echo -e "   Examples of expected values: smart-0.65, smart-0.8, smart-0.65-sub-0.5, smart-0.8-sub-0.3"
					echo -e "   Provided value: "${SEGMENTATION_ID}
					echo -e "   Segmentation condition: "${segmentation_condition}
					echo -e "   Extracted value: "${given_value}
					echo -e "   Please review column SEGMENTATION_ID in the projects master ("${MASTERPROJECTS}")."
					echo -e "   For more details on how to define segmentation see README.txt."
					echo -e ""${COL2}
					touch ./error
					exit
				fi
			}
			# check_segm_value_num "SEGM_DIFFERENCE" "${SEGM_DIFFERENCE}"
			check_segm_value_num "SEGM_MINSIZE" "${SEGM_MINSIZE}"
			check_segm_value_num "SEGM_MINSIZE_SUB" "${SEGM_MINSIZE_SUB}"
			check_segm_value_num "SEGM_MINPROBE" "${SEGM_MINPROBE}"
			check_segm_value_num "SEGM_MINPROBE_SUB" "${SEGM_MINPROBE_SUB}"
			check_segm_value_num "SEGM_MINKEEP" "${SEGM_MINKEEP}"
			check_segm_value_num "SEGM_MINLOSS" "${SEGM_MINLOSS}"
			check_segm_value_num "SEGM_MINLOSS_SUB" "${SEGM_MINLOSS_SUB}"
			check_segm_value_num "SEGM_MINLOSS_BIAL" "${SEGM_MINLOSS_BIAL}"
			check_segm_value_num "SEGM_MINGAIN" "${SEGM_MINGAIN}"
			check_segm_value_num "SEGM_MINGAIN_SUB" "${SEGM_MINGAIN_SUB}"
			check_segm_value_num "SEGM_GAP" "${SEGM_GAP}"
			check_segm_value_num "SEGM_SMOOTHPERC" "${SEGM_SMOOTHPERC}"
			check_segm_value_num "SEGM_AFDIF" "${SEGM_AFDIF}"
			check_segm_value_num "SEGM_AFSIZE" "${SEGM_AFSIZE}"
			check_segm_value_num "SEGM_AFPROBE" "${SEGM_AFPROBE}"
		fi
		#CHECKPOINT - segmentation conditions 
		#============================================================================

	fi
	###====================================================================================


	#assign databases automatically
	###====================================================================================
	if [ ${GENOME_VERSION} == "GRCh37-hg19" ]; then
		if [ `grep "chr" ${CONTIGS_SIZES} | wc -l` == 0 ]; then
			GNOMAD=${DATABASES_DIR}"GNOMAD/gnomad_2.1_GRCh37/gnomad_2.1_combine_processed_GRCh37_22XY.vcf.gz"
			CHROMOSOME_BANDS=${DATABASES_DIR}"UCSC/chr_bands/GRCh37_cytoBand.txt"
		else
			GNOMAD=${DATABASES_DIR}"GNOMAD/gnomad_2.1_hg19_liftover/gnomad_2.1_combine_processed_hg19_22XY.vcf.gz"
			CHROMOSOME_BANDS=${DATABASES_DIR}"UCSC/chr_bands/hg19_cytoBand.txt"
		fi
	fi
	if [ ${GENOME_VERSION} == "GRCh38-hg38" ]; then
		if [ `grep "chr" ${CONTIGS_SIZES} | wc -l` != 0 ]; then
	 		GNOMAD=${DATABASES_DIR}"GNOMAD/gnomad_2.1_hg38_liftover/gnomad_2.1_combine_processed_hg38_22XY.vcf.gz"
	 		CHROMOSOME_BANDS=${DATABASES_DIR}"UCSC/chr_bands/GRCh38_hg38_cytoBand.txt"
	 	else
			GNOMAD=${DATABASES_DIR}"GNOMAD/gnomad_2.1_GRCh38_liftover/gnomad_2.1_combine_processed_GRCh38_22XY.vcf.gz"
	 		CHROMOSOME_BANDS=${DATABASES_DIR}"UCSC/chr_bands/GRCh38_cytoBand.txt"
	 	fi
	fi
	if [ ${GENOME_VERSION} == "CHM13v2.0" ]; then
		GNOMAD=${DATABASES_DIR}"T2T/chm13v2.0/chm13v2.0_dbSNPv155_gnomad_processed.vcf.gz"
		CHROMOSOME_BANDS=${DATABASES_DIR}"T2T/chm13v2.0/chm13v2.0_cytoBand.txt"
	fi
	###====================================================================================



	## prepare gnomad SNP database 
	###====================================================================================

		### assign variables
		#============================================================================
		GNOMAD_ID=$(basename ${GNOMAD})
		GNOMAD_ID=`echo "${GNOMAD_ID//.vcf/}"`
		GNOMAD_ID=`echo "${GNOMAD_ID//.gz/}"`
		# GNOMAD_ID="${GNOMAD_ID%.*}"
		# GNOMAD_ID="${GNOMAD_ID%.*}"
		GNOMAD_DIR=${SUPPORT_FILES}"gnomAD/"

		GNOMAD_CONTIG_FILTER=${GNOMAD_DIR}${GNOMAD_ID}"_0_filtered_for_ref_contigs.vcf.gz"

		GNOMAD_PRESORT=${GNOMAD_DIR}${GNOMAD_ID}"_temp_presort.vcf.gz"
		GNOMAD_SORTED=${GNOMAD_DIR}${GNOMAD_ID}"_0_sorted_by_ref.vcf.gz"

		GNOMAD_CAPTURE=${GNOMAD_DIR}${GNOMAD_ID}"_1_filtered_for_capture.vcf.gz"
		GNOMAD_CAPTURE_DENOIS=${GNOMAD_DIR}${GNOMAD_ID}"_2_filter_for_capture_and_denoised_NofCTRL="${N}".vcf.gz"

		# AF_GNOMAD="0.001"
		GNOMAD_AF_TEMP_SPLIT=${GNOMAD_DIR}${GNOMAD_ID}"_1_filtered_for_af_"${AF_GNOMAD}"_temp_split.vcf"
		GNOMAD_AF_TEMP_FILTER=${GNOMAD_DIR}${GNOMAD_ID}"_1_filtered_for_af_"${AF_GNOMAD}"_temp_filter.vcf"
		GNOMAD_AF=${GNOMAD_DIR}${GNOMAD_ID}"_1_filtered_for_af_"${AF_GNOMAD}".vcf.gz"
		GNOMAD_AF_DENOIS=${GNOMAD_DIR}${GNOMAD_ID}"_2_filtered_for_af_"${AF_GNOMAD}"_and_denoised_NofCTRL="${N}".vcf.gz"

		[ ${GNOMAD_SELECTION} == "capture_filter" ] || [ ${GNOMAD_SELECTION} == "capture_filter_and_ctrl_denois" ] && GNOMAD_SELECTED_COUNT=${GNOMAD_CAPTURE} && GNOMAD_COUNT_PATTERN="_filter_capture"
		[ ${GNOMAD_SELECTION} == "af_filter" ] || [ ${GNOMAD_SELECTION} == "af_filter_and_ctrl_denois" ] && GNOMAD_SELECTED_COUNT=${GNOMAD_AF} && GNOMAD_COUNT_PATTERN="_filter_af_"${AF_GNOMAD}
		#============================================================================

		#CHECK GNOMAD
		if ! ls ${GNOMAD_SELECTED_COUNT} 1> /dev/null 2>&1; then

			### filter and sort gnomad database by reference file, if needed
			#============================================================================
			#test first if sorting needed by list of contigs
			# CONTIG_in_GNOMAD=`${GUNZIP} -c -d ${GNOMAD} | head -3000 | grep "##contig=" | awk -F"=|," '{print $3}'`
			CONTIG_in_GNOMAD=`${GUNZIP} -c -d ${GNOMAD} | awk -F "\t" '{print $1}' | uniq | grep -v "#"`
			CONTIG_in_REF=`cat ${REF_DICT} | grep "@SQ" | awk -F"\t|:" '{print $3}'`

			mkdir -p ${GNOMAD_DIR}

			#if contig lists are the same (same contigs in the same order), skip filtering for contigs and sorting of gnomAD
			if [[ "${CONTIG_in_GNOMAD}" == "${CONTIG_in_REF}" ]]; then
				GNOMAD_SORTED=${GNOMAD}
			fi

			#if contigs are not the same, filter for contigs in reference and then sort by reference
			if [[ "${CONTIG_in_GNOMAD}" != "${CONTIG_in_REF}" ]]; then

				#create list of contigs, in an order of gnomAD
				echo -e ${COL1}${BEG}${COL2}
				echo -e ${COL1}"◼︎ $(date)"${COL2}
				echo -e ${COL1}"◼︎ gnomAD SNP database - comparing contigs to reference"${COL2}
				echo -e ${COL1}${END}${COL2}
				LIST_OF_CONTIGS=""
				LIST_OF_MISSING_IN_REF=""
				for contig_in_gnomad in ${CONTIG_in_GNOMAD}; do
					if [[ ${CONTIG_in_REF[*]} =~ (^|[[:space:]])"$contig_in_gnomad"($|[[:space:]]) ]]; then
						LIST_OF_CONTIGS=`echo $LIST_OF_CONTIGS" "$contig_in_gnomad`
					else
						LIST_OF_MISSING_IN_REF=`echo $LIST_OF_MISSING_IN_REF" "$contig_in_gnomad`
					fi
				done

				#delete a space at the beggining if present (when one contig in analysis)
				LIST_OF_CONTIGS=`echo $LIST_OF_CONTIGS`
				LIST_OF_MISSING_IN_REF=`echo $LIST_OF_MISSING_IN_REF`

				#skip sorting for contigs if all contigs in gnomAD are present in reference
				if [[ "${LIST_OF_MISSING_IN_REF}" == "" ]]; then
					GNOMAD_CONTIG_FILTER=${GNOMAD}
				fi

				#filter gnomAD for contigs known in reference only
				if [[ "${LIST_OF_MISSING_IN_REF}" != "" ]]; then
					if ! ls ${GNOMAD_CONTIG_FILTER} 1> /dev/null 2>&1; then
						echo -e ${COL1}${BEG}${COL2}
						echo -e ${COL1}"◼︎ $(date)"${COL2}
						echo -e ${COL1}"◼︎ gnomAD SNP database - sorting for contigs in reference"${COL2}
						echo -e ${COL1}${END}${COL2}
						LIST_OF_CONTIGS=${LIST_OF_CONTIGS// / -L }
						${GATK} SelectVariants \
							-V ${GNOMAD} \
							-L ${LIST_OF_CONTIGS} \
							--lenient \
							-O ${GNOMAD_CONTIG_FILTER}
			 			if [[ ${CHECKPOINTS} == "yes" ]]; then 
							if ! ls ${GNOMAD_CONTIG_FILTER} 1> /dev/null 2>&1; then
								echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - gnomAD SNP PREP: file not generated"
								echo -e "   Project: "${PROJECT_ID}
								echo -e "   Stage: filtering contigs in gnomAD based on reference"
								echo -e "   File not generated: "${GNOMAD_CONTIG_FILTER}
								echo -e ""${COL2}
								touch ./error
								exit
							fi
						fi
					fi
				fi

				#sort gnomAD based on reference
				if ! ls ${GNOMAD_SORTED} 1> /dev/null 2>&1; then
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ gnomAD SNP database - sorting based on reference"${COL2}
					echo -e ${COL1}${END}${COL2}
					#decompress and delete contigs from header
					${GUNZIP} -c -d ${GNOMAD_CONTIG_FILTER} | grep -v "##contig=" | ${BGZIP} -c > ${GNOMAD_PRESORT}
					#sort by ref dict
					${PICARD} SortVcf \
						SEQUENCE_DICTIONARY=${REF_DICT} \
						I=${GNOMAD_PRESORT} \
						O=${GNOMAD_SORTED}
					rm ${GNOMAD_PRESORT}
		 			if [[ ${CHECKPOINTS} == "yes" ]]; then 
						if ! ls ${GNOMAD_SORTED} 1> /dev/null 2>&1; then
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - gnomAD SNP PREP: file not generated"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Stage: sorting SNPs in gnomAD based on reference"
							echo -e "   File not generated: "${GNOMAD_SORTED}
							echo -e ""${COL2}
							touch ./error
							exit
						fi
					fi
				fi
			fi
			#============================================================================

			### filter gnomad database by capture, if needed
			#============================================================================
			if [ ${GNOMAD_SELECTION} == "capture_filter" ] || [ ${GNOMAD_SELECTION} == "capture_filter_and_ctrl_denois" ]; then 
				if ! ls ${GNOMAD_CAPTURE} 1> /dev/null 2>&1; then
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ gnomAD SNP database - filtering for capture - "${CAPTURE_ID}${COL2}
					echo -e ${COL1}${END}${COL2}
					mkdir -p ${GNOMAD_DIR}
					${GATK} SelectVariants \
						-R ${REF} \
						-V ${GNOMAD_SORTED} \
						-L ${INTERVALS_PROCESSED} \
						--lenient \
						-O ${GNOMAD_CAPTURE}
		 			if [[ ${CHECKPOINTS} == "yes" ]]; then 
						if ! ls ${GNOMAD_CAPTURE} 1> /dev/null 2>&1; then
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - gnomAD SNP PREP: file not generated"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Stage: filtering gnomAD SNPs for capture"
							echo -e "   File not generated: "${GNOMAD_CAPTURE}
							echo -e ""${COL2}
							touch ./error
							exit
						fi
					fi
				fi
			fi
			#============================================================================

			### filter gnomad database by AF > spec freq, if needed
			#============================================================================
			if [ ${GNOMAD_SELECTION} == "af_filter" ] || [ ${GNOMAD_SELECTION} == "af_filter_and_ctrl_denois" ]; then 
				if ! ls ${GNOMAD_AF} 1> /dev/null 2>&1; then
					#split multiallelic sites and indes - 25min
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ filtering gnomad database for AF > "${AF_GNOMAD}" - splitting multiallelic sites"${COL2}
					echo -e ${COL1}${END}${COL2}
					mkdir -p ${GNOMAD_DIR}
					${BCFTOOLS} norm -m -both --do-not-normalize --check-ref ws -f ${REF} \
						-Ov -o ${GNOMAD_AF_TEMP_SPLIT} \
						${GNOMAD_SORTED}
					${GATK} IndexFeatureFile \
	 					-I ${GNOMAD_AF_TEMP_SPLIT}
					#filter for AF - 8min
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ filtering gnomad database for AF > "${AF_GNOMAD}" - filtering for AF"${COL2}
					echo -e ${COL1}${END}${COL2}
					${GATK} SelectVariants \
						-R ${REF} \
						-V ${GNOMAD_AF_TEMP_SPLIT} \
						-L ${INTERVALS_PROCESSED} \
						-select "AF > "${AF_GNOMAD} \
						--lenient \
						-O ${GNOMAD_AF_TEMP_FILTER}
					#merge multiallelic back and indes - 3min
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ filtering gnomad database for AF > "${AF_GNOMAD}" - merging multiallelic sites back"${COL2}
					echo -e ${COL1}${END}${COL2}
					${BCFTOOLS} norm -m +any --do-not-normalize --check-ref ws -f ${REF} \
						-Oz -o ${GNOMAD_AF} \
						${GNOMAD_AF_TEMP_FILTER}
					${GATK} IndexFeatureFile \
	 					-I ${GNOMAD_AF}
					rm ${GNOMAD_AF_TEMP_SPLIT}*
					rm ${GNOMAD_AF_TEMP_FILTER}*
		 			if [[ ${CHECKPOINTS} == "yes" ]]; then 
						if ! ls ${GNOMAD_AF} 1> /dev/null 2>&1; then
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - gnomAD SNP PREP: file not generated"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Stage: filtering gnomAD SNPs for AF > "${AF_GNOMAD}
							echo -e "   File not generated: "${GNOMAD_AF}
							echo -e ""${COL2}
							touch ./error
							exit
						fi
					fi
				fi
			fi
			#============================================================================

		#CHECK GNOMAD
		fi

	###====================================================================================


	### delete controls that were exclueded (yes -> no)
	###====================================================================================

		#CHECK - CTRL YES 
		if [ ${CTRL} == "yes" ]; then 
			cat ${MASTERCONTROLS} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" '{if ($1=="no" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION)print}' | while read -r line || [[ -n "$line" ]]; do
			 	MAIN_ID=`extract_from_master 'MAIN_ID' ${MASTERCONTROLS}`
				CTRL_ID=`extract_from_master 'CTRL_ID' ${MASTERCONTROLS}`
				CTRL_SEX=`extract_from_master 'CTRL_SEX' ${MASTERCONTROLS}`
				RESULTS_DATA_PATTERN=${MAIN_ID}"_"${CTRL_ID}"_"${CTRL_SEX}
				if ls ${CTRL_COUNT_HDF_DIR}${RESULTS_DATA_PATTERN}"_counts.hdf5" 1> /dev/null 2>&1; then
					rm ${CTRL_COUNT_HDF_DIR}${RESULTS_DATA_PATTERN}"_counts.hdf5"
				fi
				if ls ${CTRL_COUNT_TSV_DIR}${RESULTS_DATA_PATTERN}"_counts.tsv" 1> /dev/null 2>&1; then
					rm ${CTRL_COUNT_TSV_DIR}${RESULTS_DATA_PATTERN}"_counts.tsv"
				fi
				if ls ${CTRL_MAF_DIR}${RESULTS_DATA_PATTERN}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".rds" 1> /dev/null 2>&1; then
					rm ${CTRL_MAF_DIR}${RESULTS_DATA_PATTERN}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".rds"
				fi
				if ls ${CTRL_SEX_DIR}${RESULTS_DATA_PATTERN}"_sextest"* 1> /dev/null 2>&1; then
					rm ${CTRL_SEX_DIR}${RESULTS_DATA_PATTERN}"_sextest"*
				fi
				if ls ${CTRL_SEX_DIR}"results"* 1> /dev/null 2>&1; then
					rm ${CTRL_SEX_DIR}"results"*
				fi
			done
		#CHECK - CTRL YES 
		fi
		#============================================================================
	
	###====================================================================================



	### sex test - if any control does not give expected sex, exit
	###====================================================================================

		#CHECK - CTRL YES 
		if [ ${CTRL} == "yes" ]; then

			### skip if sex test is not required
			#============================================================================
			if [ ${CTRL_SEX_TEST} == "no" ]; then
				echo -e ${COL1}${BEG}${COL2}
				echo -e ${COL1}"◼︎ $(date)"${COL2}
				echo -e ${COL1}"◼︎ skipping sex test of controls"${COL2}
				echo -e ${COL1}${END}${COL2}
			fi
			#============================================================================

			### loop the controls, collect read count, allelic count and text sex
			#============================================================================
			if [ ${CTRL_SEX_TEST} == "default" ] || [ ${CTRL_SEX_TEST} == "custom" ]; then
			
				#default setting, futher selection based on genome version and if there is "chr" in GRCh37-hg19
				[ ${CTRL_SEX_TEST} == "default" ] && [ ${GENOME_VERSION} == "GRCh37-hg19" ] && [ `grep "chr" ${CONTIGS_SIZES} | wc -l` == 0 ] && SEX_TEST_SETTING="DEFAULT_1"
				[ ${CTRL_SEX_TEST} == "default" ] && [ ${GENOME_VERSION} == "GRCh37-hg19" ] && [ `grep "chr" ${CONTIGS_SIZES} | wc -l` != 0 ] && SEX_TEST_SETTING="DEFAULT_2"
				[ ${CTRL_SEX_TEST} == "default" ] && [ ${GENOME_VERSION} == "GRCh38-hg38" ] && SEX_TEST_SETTING="DEFAULT_3"
				[ ${CTRL_SEX_TEST} == "default" ] && [ ${GENOME_VERSION} == "CHM13v2.0" ] && SEX_TEST_SETTING="DEFAULT_4"

				[ ${CTRL_SEX_TEST} == "custom" ] && SEX_TEST_SETTING=${PROJECT_ID}

				MASTERSEX_FILTER=${MASTERSEX}"_temp.bed"
				cat ${MASTERSEX} | awk -F"\t" -v SEX_TEST_SETTING="${SEX_TEST_SETTING}" '{if ($5==SEX_TEST_SETTING)print $1"\t"$2"\t"$3"\t"$4}' > ${MASTERSEX_FILTER}

				#CHECKPOINT - sex master
				#============================================================================
				if [[ ${CHECKPOINTS} == "yes" ]]; then 
					#test if filtered master sex is not empty
					NofLINES=`cat ${MASTERSEX_FILTER} | wc -l`
					if [[ $NofLINES == 0 ]]; then
						echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - SEX TEST: could not find any setting defined in the sex master ("${MASTERSAMPLES}") to perform sex test"
						echo -e "   Project: "${PROJECT_ID}
						echo -e "   Pattern that was looked for in the sex master file column PROJECT_ID: "${SEX_TEST_SETTING}
						echo -e "   If you are trying to use custom setting for the sex test, make sure that the PROJECT_ID from the projects master is same as the PROJECT_ID within the sex master."
						echo -e "   To unable sex test to be performed, change column CTRL_SEX_TEST in the projects master ("${MASTERPROJECTS}") to no."
						echo -e ""${COL2}
						touch ./error
						exit
					fi
					#test that one line was recognized for both control and test gene
					TEST_GENE_N=`cat ${MASTERSEX_FILTER} | awk -F"\t" '{if ($4=="TEST_GENE")print $1":"$2"-"$3}' | wc -l`
					CTRL_GENE_N=`cat ${MASTERSEX_FILTER} | awk -F"\t" '{if ($4=="CTRL_GENE")print $1":"$2"-"$3}' | wc -l`

					if [[ $TEST_GENE_N == 1 ]] && [[ $CTRL_GENE_N == 1 ]]; then true; else
						echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - SEX TEST: found unexpected number of lines for tested and/or control gene defined in sex master ("${MASTERSAMPLES}")"
						echo -e "   Project: "${PROJECT_ID}
						echo -e "   Tested gene (TEST_GENE in column GENE_TYPE) expected to have 1 line and has "${TEST_GENE_N}.
						echo -e "   Control gene (CTRL_GENE in column GENE_TYPE) expected to have 1 line and has "${CTRL_GENE_N}.
						echo -e "   To unable sex test to be performed, change column CTRL_SEX_TEST in the projects master ("${MASTERPROJECTS}") to no."
						echo -e ""${COL2}
						touch ./error
						exit
					fi
					#test if contigs is present in reference file and if tested gene is on chromosome Y
					TEST_GENE_CONTIG=`cat ${MASTERSEX_FILTER} | awk -F"\t" '{if ($4=="TEST_GENE")print $1}'`
					CTRL_GENE_CONTIG=`cat ${MASTERSEX_FILTER} | awk -F"\t" '{if ($4=="CTRL_GENE")print $1}'`
					CONTIG_in_REF=`cat ${REF_DICT} | grep "@SQ" | awk -F"\t|:" '{print $3}'`
					if [[ ${CONTIG_in_REF} =~ (^|[[:space:]])"${TEST_GENE_CONTIG}"($|[[:space:]]) ]] && [[ ${CONTIG_in_REF} =~ (^|[[:space:]])"${CTRL_GENE_CONTIG}"($|[[:space:]]) ]]; then true; else
						echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - SEX TEST: contig for tested and/or control gene is not present in the reference file"
						echo -e "   Project: "${PROJECT_ID}
						echo -e "   Given contig for tested gene: "${TEST_GENE_CONTIG}
						echo -e "   Given contig for control gene: "${CTRL_GENE_CONTIG}
						echo -e "   Contigs in reference file: "${CONTIG_in_REF}
						echo -e "   Please correct this in the sex master ("${MASTERSAMPLES}")."
						echo -e "   To unable sex test to be performed, change column CTRL_SEX_TEST in the projects master ("${MASTERPROJECTS}") to no."
						echo -e ""${COL2}
						touch ./error
						exit
					fi
					if [[ $TEST_GENE_CONTIG == "Y" ]] || [[ $TEST_GENE_CONTIG == "chrY" ]]; then true; else
						echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - SEX TEST: contig for tested gene is not on chromosome Y"
						echo -e "   Project: "${PROJECT_ID}
						echo -e "   Given contig for tested gene: "${TEST_GENE_CONTIG}
						echo -e "   This test requires chromosome Y locus to be tested. If no locus on chromosome Y is captured, sex test cannot be performed."
						echo -e "   To unable sex test to be performed, change column CTRL_SEX_TEST in the projects master ("${MASTERPROJECTS}") to no."
						echo -e ""${COL2}
						touch ./error
						exit
					fi
				fi
				#CHECKPOINT - sex master
				#============================================================================


				#collect coverage data
			 	cat ${MASTERCONTROLS} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION)print}' | while read -r line || [[ -n "$line" ]]; do
					
					#prevent looping if error occured during previous round
					if [[ ${CHECKPOINTS} == "yes" ]]; then 
						if ls ./error 1> /dev/null 2>&1; then exit; fi
					fi

					MAIN_ID=`extract_from_master 'MAIN_ID' ${MASTERCONTROLS}`
					CTRL_ID=`extract_from_master 'CTRL_ID' ${MASTERCONTROLS}`
					CTRL_SEX=`extract_from_master 'CTRL_SEX' ${MASTERCONTROLS}`
					CTRL_BAM_PATH=`extract_from_master 'CTRL_PATH_TO_BAM' ${MASTERCONTROLS}`

					[ ${WAY_TO_BAM} == "find_in_dir" ] && BAMFILE=`find ${CTRL_BAM_DIR} -type f -name ${CTRL_ID}*${CTRL_BAM_PATTERN}`
					[ ${WAY_TO_BAM} == "absolute" ] && BAMFILE=${CTRL_BAM_DIR}${CTRL_BAM_PATH}

					RESULTS_DATA_PATTERN=${MAIN_ID}"_"${CTRL_ID}"_"${CTRL_SEX}
					ECHO_PATTERN=${MAIN_ID}" - "${CTRL_ID}

					### sex test - prepare read counts based on SRY and GAPDH (or alternative genes) coverage
					OUTPUT=${CTRL_SEX_DIR}${RESULTS_DATA_PATTERN}"_sextest.cov"
					if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
						echo -e ${COL1}${BEG}${COL2}
						echo -e ${COL1}"◼︎ $(date)"${COL2}
						echo -e ${COL1}"◼︎ controls prep - "${ECHO_PATTERN}": preparing read counts for sex test"${COL2}
						mkdir -p ${CTRL_SEX_DIR}
						#extract positions of interest to bam and index bam
						echo -e ${COL1}"  ... creating subset bam"${COL2}
						## this is much much much slower:
						# ${SAMTOOLS} view -b -L ${MASTERSEX_FILTER} ${BAMFILE} > ${CTRL_SEX_DIR}${RESULTS_DATA_PATTERN}"_sextest.bam"
						## this is quicker ${SAMTOOLS} view option, but gives warnings with older indexes, so just ignore it:
						TEST_GENE=`cat ${MASTERSEX_FILTER} | awk -F"\t" '{if ($4=="TEST_GENE")print $1":"$2"-"$3}'`
						CTRL_GENE=`cat ${MASTERSEX_FILTER} | awk -F"\t" '{if ($4=="CTRL_GENE")print $1":"$2"-"$3}'`
						${SAMTOOLS} view -b ${BAMFILE} ${TEST_GENE} ${CTRL_GENE} > ${CTRL_SEX_DIR}${RESULTS_DATA_PATTERN}"_sextest.bam"
						${SAMTOOLS} index ${CTRL_SEX_DIR}${RESULTS_DATA_PATTERN}"_sextest.bam"
						#get coverage for regions of interest
						echo -e ${COL1}"  ... and processing coverage"${COL2}
						${BEDTOOLS} coverage -a ${MASTERSEX_FILTER} -b ${CTRL_SEX_DIR}${RESULTS_DATA_PATTERN}"_sextest.bam" > ${OUTPUT}
						echo -e ${COL1}${END}${COL2}
						if [[ ${CHECKPOINTS} == "yes" ]]; then 
							if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
								echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
								echo -e "   Project: "${PROJECT_ID}
								echo -e "   Main ID: "${MAIN_ID}
								echo -e "   Control ID: "${CTRL_ID}
								echo -e "   Stage: sex test"
								echo -e "   File not generated: "${OUTPUT}
								echo -e ""${COL2}
								touch ./error
								exit
							fi
						fi
					fi
				done

				rm ${MASTERSEX_FILTER}
				
				#prevent further processing if error occured during previous step
				if [[ ${CHECKPOINTS} == "yes" ]]; then 
					if ls ./error 1> /dev/null 2>&1; then exit; fi
				fi

				### process all cov files 
				echo -e ${COL1}${BEG}${COL2}
				echo -e ${COL1}"◼︎ $(date)"${COL2}
				echo -e ${COL1}"◼︎ controls prep - testing sex"${COL2}
				if ls ${CTRL_SEX_DIR}"results"* 1> /dev/null 2>&1; then
					rm ${CTRL_SEX_DIR}"results"*
				fi
				Rscript --vanilla ${R}"ctrl_sex_test.R" \
						${CTRL_SEX_DIR}
				### check size of the error file and if it contains anything, exit the running
				FILESIZE=$(wc -c < ${CTRL_SEX_DIR}"results_wrong.tsv")
				if [ ${FILESIZE} = 0 ]; then
					echo -e ${COL1}"◼︎ ... controls have expected sex"${COL2}
					echo -e ${COL1}${END}${COL2}
				fi
				if [ ${FILESIZE} != 0 ]; then
					echo ""
					echo -e ${COL1_WARN}"◼︎ !!! CHECK SEX OF CONTROLS !!!"${COL2}
					echo -e ${COL1_WARN}
					cat ${CTRL_SEX_DIR}"results_wrong.tsv"
					echo -e ${COL2}
					echo -e ${COL1}${END}${COL2}
					touch ./error
					exit
				fi

			fi
			#============================================================================
		
		#CHECK - CTRL YES 
		fi
	###====================================================================================


	### process controls: read count, allelic count, PoN, denoised, segmentation and QC
	###====================================================================================

		#CHECK - CTRL YES 
		if [ ${CTRL} == "yes" ]; then

			### loop the controls, collect read count and allelic count
			#============================================================================
			#LOOP - CTRLS
			cat ${MASTERCONTROLS} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION)print}' | while read -r line || [[ -n "$line" ]]; do

				#prevent looping if error occured during previous round
				if [[ ${CHECKPOINTS} == "yes" ]]; then 
					if ls ./error 1> /dev/null 2>&1; then exit; fi
				fi

				MAIN_ID=`extract_from_master 'MAIN_ID' ${MASTERCONTROLS}`
				CTRL_ID=`extract_from_master 'CTRL_ID' ${MASTERCONTROLS}`
				CTRL_SEX=`extract_from_master 'CTRL_SEX' ${MASTERCONTROLS}`
				CTRL_BAM_PATH=`extract_from_master 'CTRL_PATH_TO_BAM' ${MASTERCONTROLS}`

				[ ${WAY_TO_BAM} == "find_in_dir" ] && BAMFILE=`find ${CTRL_BAM_DIR} -type f -name ${CTRL_ID}*${CTRL_BAM_PATTERN}`
				[ ${WAY_TO_BAM} == "absolute" ] && BAMFILE=${CTRL_BAM_DIR}${CTRL_BAM_PATH}

				RESULTS_DATA_PATTERN=${MAIN_ID}"_"${CTRL_ID}"_"${CTRL_SEX}
				ECHO_PATTERN=${MAIN_ID}" - "${CTRL_ID}

				### collect read count for controls, for PoN and for CN freq
				#LOOP - READ COUNTS FORMAT
				for FILE_TYPE in "HDF5" "TSV"; do
					[ ${FILE_TYPE} == "HDF5" ] && FILE_TYPE_END="hdf5"
					[ ${FILE_TYPE} == "TSV" ] && FILE_TYPE_END="tsv"
					[ ${FILE_TYPE} == "HDF5" ] && CTRL_COUNT_DIR=${CTRL_COUNT_HDF_DIR}
					[ ${FILE_TYPE} == "TSV" ] && CTRL_COUNT_DIR=${CTRL_COUNT_TSV_DIR}
					OUTPUT=${CTRL_COUNT_DIR}${RESULTS_DATA_PATTERN}"_counts."${FILE_TYPE_END}
					if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
						echo -e ${COL1}${BEG}${COL2}
						echo -e ${COL1}"◼︎ $(date)"${COL2}
						echo -e ${COL1}"◼︎ controls prep - "${ECHO_PATTERN}": collecting read count - "${FILE_TYPE_END}${COL2}
						echo -e ${COL1}${END}${COL2}
						mkdir -p ${CTRL_COUNT_DIR}
						if [ ${READ_COUNT_FILTER} == "F" ]; then
							${GATK} CollectReadCounts \
								-R ${REF} \
								-I ${BAMFILE} \
								-L ${INTERVALS_PROCESSED} \
								--interval-merging-rule OVERLAPPING_ONLY \
								--format ${FILE_TYPE} \
								-O ${OUTPUT}
						fi
						if [ ${READ_COUNT_FILTER} == "DF" ]; then
							${GATK} CollectReadCounts \
								-R ${REF} \
								-DF "MappingQualityReadFilter" \
								-I ${BAMFILE} \
								-L ${INTERVALS_PROCESSED} \
								--interval-merging-rule OVERLAPPING_ONLY \
								--format ${FILE_TYPE} \
								-O ${OUTPUT}
						fi
						if [[ ${CHECKPOINTS} == "yes" ]]; then 
							if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
								echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
								echo -e "   Project: "${PROJECT_ID}
								echo -e "   Main ID: "${MAIN_ID}
								echo -e "   Control ID: "${CTRL_ID}
								echo -e "   Stage: collection of read counts"
								echo -e "   File not generated: "${OUTPUT}
								echo -e ""${COL2}
								touch ./error
								exit
							fi
						fi
					fi
				done

				## collect allelic counts for denoising SNP profile and convert to rds (-RF GoodCigarReadFilter \)
				OUTPUT=${CTRL_MAF_DIR}${RESULTS_DATA_PATTERN}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".rds"
				if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ controls prep - "${ECHO_PATTERN}": collecting allelic count"${COL2}
					echo -e ${COL1}${END}${COL2}
					mkdir -p ${CTRL_MAF_DIR}
					${GATK} CollectAllelicCounts \
						-R ${REF} \
						-I ${BAMFILE} \
				 		-L ${GNOMAD_SELECTED_COUNT} \
				 		-O ${OUTPUT}".tsv"
		 			echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ controls prep - "${ECHO_PATTERN}": converting allelic counts file to rds"${COL2} 
					echo -e ${COL1}${END}${COL2}
					Rscript --vanilla ${R}convert_alleliccounts_to_rds.R \
							${OUTPUT}".tsv" \
							${OUTPUT}
					rm ${OUTPUT}".tsv"
		 			if [[ ${CHECKPOINTS} == "yes" ]]; then 
						if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Main ID: "${MAIN_ID}
							echo -e "   Control ID: "${CTRL_ID}
							echo -e "   Stage: collection of allelic counts"
							echo -e "   File not generated: "${OUTPUT}
							echo -e ""${COL2}
							touch ./error
							exit
						fi
					fi
			 	fi

			#LOOP - CTRLS
			done

			#prevent further processing if error occured during previous step
			if [[ ${CHECKPOINTS} == "yes" ]]; then 
				if ls ./error 1> /dev/null 2>&1; then exit; fi
			fi
			#============================================================================

			### create PoN for each control with appropriate sex and denois each control
			#============================================================================
			#LOOP - CTRLS
			cat ${MASTERCONTROLS} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION)print}' | while read -r line || [[ -n "$line" ]]; do

				#prevent looping if error occured during previous round
				if [[ ${CHECKPOINTS} == "yes" ]]; then 
					if ls ./error 1> /dev/null 2>&1; then exit; fi
				fi

				MAIN_ID=`extract_from_master 'MAIN_ID' ${MASTERCONTROLS}`
				CTRL_ID=`extract_from_master 'CTRL_ID' ${MASTERCONTROLS}`
				CTRL_SEX=`extract_from_master 'CTRL_SEX' ${MASTERCONTROLS}`

				RESULTS_DATA_PATTERN=${MAIN_ID}"_"${CTRL_ID}"_"${CTRL_SEX}
				ECHO_PATTERN=${MAIN_ID}" - "${CTRL_ID}
			
				[ ${CTRL_PON_SEX_SELECT} == "M" ] && PON_SEX_PATTERN="M" && PON_SEX_ID="M"
				[ ${CTRL_PON_SEX_SELECT} == "F" ] && PON_SEX_PATTERN="F" && PON_SEX_ID="F"
				[ ${CTRL_PON_SEX_SELECT} == "mixed" ] && PON_SEX_PATTERN="*" && PON_SEX_ID="mixed"
				[ ${CTRL_PON_SEX_SELECT} == "matched" ] && PON_SEX_PATTERN=${CTRL_SEX} && PON_SEX_ID=${CTRL_SEX}

				#generate list of hdf5 controls from controls master (list, add -I and delet newline \n)
				if [ ${CTRL_PON_SEX_SELECT} == "mixed" ]; then
					cat ${MASTERCONTROLS} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" -v MAIN_ID="${MAIN_ID}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION && $5!=MAIN_ID)print}' > ${MASTERCONTROLS}"_temp"
				fi
				if [ ${CTRL_PON_SEX_SELECT} == "M" ] || [ ${CTRL_PON_SEX_SELECT} == "F" ] || [ ${CTRL_PON_SEX_SELECT} == "matched" ]; then
					cat ${MASTERCONTROLS} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" -v MAIN_ID="${MAIN_ID}" -v PON_SEX_PATTERN="${PON_SEX_PATTERN}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION && $5!=MAIN_ID && $7==PON_SEX_PATTERN)print}' > ${MASTERCONTROLS}"_temp"
				fi
				hdflist=`awk -F'\t' -v CTRL_COUNT_HDF_DIR="${CTRL_COUNT_HDF_DIR}" '{print CTRL_COUNT_HDF_DIR $5 "_" $6 "_" $7 "_counts.hdf5"}' ${MASTERCONTROLS}"_temp" | sed 's/^/-I /'`
				NRealx=`cat ${MASTERCONTROLS}"_temp" | wc -l`
				NReal=$(echo ${NRealx} | tr -d ' ')
				if ls ${MASTERCONTROLS}"_temp" 1> /dev/null 2>&1; then
					rm ${MASTERCONTROLS}"_temp"
				fi

				#prepare PoN
				if [ ${NReal} != 0 ]; then
					OUTPUT=${CTRL_PoN_DIR}${MAIN_ID}"_PoN-"${PON_SEX_ID}".hdf5"
					if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
						echo -e ${COL1}${BEG}${COL2}
						echo -e ${COL1}"◼︎ $(date)"${COL2}
						echo -e ${COL1}"◼︎ controls prep - "${ECHO_PATTERN}": preparing PoN-"${PON_SEX_ID}${COL2}
						echo -e ${COL1}${END}${COL2}
						mkdir -p ${CTRL_PoN_DIR}
						echo ${hdflist} > ${CTRL_PoN_DIR}${MAIN_ID}"_PoN-"${PON_SEX_ID}"_list.txt"
						${GATK} CreateReadCountPanelOfNormals \
							${hdflist} \
							--annotated-intervals ${INTERVALS_ANNOTATED} \
							-O ${OUTPUT}
			 			if [[ ${CHECKPOINTS} == "yes" ]]; then 
							if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
								echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
								echo -e "   Project: "${PROJECT_ID}
								echo -e "   Main ID: "${MAIN_ID}
								echo -e "   Control ID: "${CTRL_ID}
								echo -e "   Stage: panel of normals"
								echo -e "   File not generated: "${OUTPUT}
								echo -e ""${COL2}
								touch ./error
								exit
							fi
						fi
					fi
				fi
				
				#denoise CN data
				OUTPUT=${CTRL_DENOIS_DIR}${RESULTS_DATA_PATTERN}"_denoisedCR_PoN-"${PON_SEX_ID}".tsv"
				if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
					#denois CN data when controls available
					if [ ${NReal} != 0 ]; then
						echo -e ${COL1}${BEG}${COL2}
						echo -e ${COL1}"◼︎ $(date)"${COL2}
						echo -e ${COL1}"◼︎ controls prep - "${ECHO_PATTERN}": preparing denoising CN data"${COL2}
						echo -e ${COL1}${END}${COL2}
						mkdir -p ${CTRL_STANDARD_DIR}
						mkdir -p ${CTRL_DENOIS_DIR}
						${GATK} DenoiseReadCounts \
							-I ${CTRL_COUNT_HDF_DIR}${RESULTS_DATA_PATTERN}"_counts.hdf5" \
							--count-panel-of-normals ${CTRL_PoN_DIR}${MAIN_ID}"_PoN-"${PON_SEX_ID}".hdf5" \
							--standardized-copy-ratios ${CTRL_STANDARD_DIR}${RESULTS_DATA_PATTERN}"_standardizedCR_PoN-"${PON_SEX_ID}".tsv" \
							--denoised-copy-ratios ${OUTPUT}
					fi
					#denoise CN data by GC correction only, when no control available
					if [ ${NReal} == 0 ]; then
						echo -e ${COL1}${BEG}${COL2}
						echo -e ${COL1}"◼︎ $(date)"${COL2}
						echo -e ${COL1}"◼︎ controls prep - "${ECHO_PATTERN}": preparing denoising CN data"${COL2}
						echo -e ${COL1}${END}${COL2}
						mkdir -p ${CTRL_STANDARD_DIR}
						mkdir -p ${CTRL_DENOIS_DIR}
						${GATK} DenoiseReadCounts \
							-I ${CTRL_COUNT_HDF_DIR}${RESULTS_DATA_PATTERN}"_counts.hdf5" \
							--annotated-intervals ${INTERVALS_ANNOTATED} \
							--standardized-copy-ratios ${CTRL_STANDARD_DIR}${RESULTS_DATA_PATTERN}"_standardizedCR_PoN-"${PON_SEX_ID}".tsv" \
							--denoised-copy-ratios ${OUTPUT}
					fi
					if [[ ${CHECKPOINTS} == "yes" ]]; then 
						if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Main ID: "${MAIN_ID}
							echo -e "   Control ID: "${CTRL_ID}
							echo -e "   Stage: CN standardization and denoising"
							echo -e "   File not generated: "${OUTPUT}
							echo -e ""${COL2}
							touch ./error
							exit
						fi
					fi
				fi

				#segmentation
				if [[ ${QUICK_MODE} != "Q" ]]; then
					OUTPUT=${CTRL_SEGMENT_DIR}${RESULTS_DATA_PATTERN}"_SEGMENTS_PoN-"${PON_SEX_ID}"_"${SEGMENTATIONCONDITIONS}"_segmentation.tsv"
					if ! ls ${OUTPUT} 1> /dev/null 2>&1; then 
						echo -e ${COL1}${BEG}${COL2}
						echo -e ${COL1}"◼︎ $(date)"${COL2}
						echo -e ${COL1}"◼︎ controls prep - "${ECHO_PATTERN}": segmentation and quality control"${COL2} 
						echo -e ${COL1}${END}${COL2}
						mkdir -p ${CTRL_SEGMENT_DIR}
						mkdir -p ${CTRL_QC_DIR}	
						
						#generate segmentation condition if smart way is required
						if [[ ${SEGMENTATION_ID} == "smart-"* ]]; then	
							SMART_SEGMENTATION=`Rscript --vanilla ${R}"segmentation_conditions_smart.R" \
							${SEGMENTATION_ID} \
							${CTRL_DENOIS_DIR}${RESULTS_DATA_PATTERN}"_denoisedCR_PoN-"${PON_SEX_ID}".tsv" \
							${INTERVALS_ANNOTATED} \
							${SEQ_TYPE}`

							#get rid off the R format
							SMART_SEGMENTATION=$(echo ${SMART_SEGMENTATION} | awk '{print $2}')
							SMART_SEGMENTATION=$(echo ${SMART_SEGMENTATION} | tr -d '\"')
							#extract values that depends on denoised data
							SEGM_DIFFERENCE=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $1}') 
		
							SEGMENTATION_FULL=$SEGM_DIFFERENCE$x$SEGM_MINSIZE$x$SEGM_MINSIZE_SUB$x$SEGM_MINPROBE$x$SEGM_MINPROBE_SUB$x$SEGM_MINKEEP$x$SEGM_MINLOSS$x$SEGM_MINLOSS_SUB$x$SEGM_MINLOSS_BIAL$x$SEGM_MINGAIN$x$SEGM_MINGAIN_SUB$x$SEGM_GAP$x$SEGM_SMOOTHPERC$x$SEGM_AFDIF$x$SEGM_AFSIZE$x$SEGM_AFPROBE
						
							#CHECKPOINT - segmentation conditions 
							#============================================================================
							if [[ ${CHECKPOINTS} == "yes" ]]; then
								#function to check if values values are numerical
								function check_segm_value_num {
									segmentation_condition=$1
									given_value=$2
									numbers='^[+-]?[0-9]+([.][0-9]+)?$'
									if [[ $given_value =~ $numbers ]]; then true; else 
										echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - SEGMENTATION CONDITIONS: some condition with unexpected value"
										echo -e "   Project: "${PROJECT_ID}
										echo -e "   Segmentation conditions were indicated to be defined under smart segmentation, but format likely not provided correctly."
										echo -e "   Examples of expected values: smart-0.65, smart-0.8, smart-0.65-sub-0.5, smart-0.8-sub-0.3"
										echo -e "   Provided value: "${SEGMENTATION_ID}
										echo -e "   Segmentation condition: "${segmentation_condition}
										echo -e "   Extracted value: "${given_value}
										echo -e "   Please review column SEGMENTATION_ID in the projects master ("${MASTERPROJECTS}")."
										echo -e "   For more details on how to define segmentation see README.txt."
										echo -e ""${COL2}
										touch ./error
										exit
									fi
								}
								check_segm_value_num "SEGM_DIFFERENCE" "${SEGM_DIFFERENCE}"
							fi
							#CHECKPOINT - segmentation conditions 
							#============================================================================
						fi

						#segmentation
						Rscript --vanilla ${R}"segmentation_and_qc.R" \
							${CHROMOSOME_BANDS} \
							${CAPTURE_ID} \
							${PROJECT_ID} \
							${PROJECT_TYPE} \
							${CTRL_COUNT_TSV_DIR}${RESULTS_DATA_PATTERN}"_counts.tsv" \
							${CTRL_MAF_DIR}${RESULTS_DATA_PATTERN}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".rds" \
							${CTRL_STANDARD_DIR}${RESULTS_DATA_PATTERN}"_standardizedCR_PoN-"${PON_SEX_ID}".tsv" \
							${CTRL_DENOIS_DIR}${RESULTS_DATA_PATTERN}"_denoisedCR_PoN-"${PON_SEX_ID}".tsv" \
							${MAIN_ID} \
							${CTRL_ID} \
							"control" \
							${CTRL_SEX} \
							${PON_SEX_ID} \
							${N} \
							${NReal} \
							${SEGM_DIFFERENCE} \
							${SEGM_MINSIZE} \
							${SEGM_MINSIZE_SUB} \
							${SEGM_MINPROBE} \
							${SEGM_MINPROBE_SUB} \
							${SEGM_MINKEEP} \
							${SEGM_MINLOSS} \
							${SEGM_MINLOSS_SUB} \
							${SEGM_MINLOSS_BIAL} \
							${SEGM_MINGAIN} \
							${SEGM_MINGAIN_SUB} \
							${SEGM_GAP} \
							${SEGM_SMOOTHPERC} \
							${SEGM_AFDIF} \
							${SEGM_AFSIZE} \
							${SEGM_AFPROBE} \
							${SEGMENTATION_ID}"_"${SEGMENTATION_FULL} \
							${CTRL_QC_DIR}"NofCTRL="${N}"_qc_summary.tsv" \
							${OUTPUT}
					
						if [[ ${CHECKPOINTS} == "yes" ]]; then 
							if ls "status_ok" 1> /dev/null 2>&1; then
								rm "status_ok"
							else
								echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: R script was halted"
								echo -e "   Project: "${PROJECT_ID}
								echo -e "   Main ID: "${MAIN_ID}
								echo -e "   Control ID: "${CTRL_ID}
								echo -e "   Stage: data segmentation"
								echo -e "   File not generated: "${OUTPUT}
								echo -e "   R script: segmentation_and_qc.R"
								echo -e ""${COL2}
								touch ./error
								exit
							fi
							if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
								echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
								echo -e "   Project: "${PROJECT_ID}
								echo -e "   Main ID: "${MAIN_ID}
								echo -e "   Control ID: "${CTRL_ID}
								echo -e "   Stage: data segmentation"
								echo -e "   File not generated: "${OUTPUT}
								echo -e ""${COL2}
								touch ./error
								exit
							fi
						fi
					
					fi
				fi
					
			#LOOP - CTRLS
			done

			#prevent further processing if error occured during previous step
			if [[ ${CHECKPOINTS} == "yes" ]]; then 
				if ls ./error 1> /dev/null 2>&1; then exit; fi
			fi
			#============================================================================

		#CHECK - CTRL YES 
		fi
	###====================================================================================

	### use control denoised CN to prepare pre-file and process by segmentation conditions
	###====================================================================================
	
		#CHECK - CTRL YES 
		if [ ${CTRL} == "yes" ]; then

			#CHECK - CTRL CN REQUIRED
			if [ ${CN_FREQ_IN_CTRLS} == "yes" ]; then

				[ ${CTRL_PON_SEX_SELECT} == "M" ] && PON_SEX_PATTERN="_denoisedCR_PoN-M"
				[ ${CTRL_PON_SEX_SELECT} == "F" ] && PON_SEX_PATTERN="_denoisedCR_PoN-F"
				[ ${CTRL_PON_SEX_SELECT} == "mixed" ] && PON_SEX_PATTERN="_denoisedCR_PoN-mixed"
				[ ${CTRL_PON_SEX_SELECT} == "matched" ] && PON_SEX_PATTERN="M_denoisedCR_PoN-M|F_denoisedCR_PoN-F"

				CTRL_CN_FILE=${CTRL_CN_DIR}"ctrl_denois_CN_NofCTRL="${N}"_PoN-"${CTRL_PON_SEX_SELECT}".rds"
				CTRL_CN_FILE_SEGMENT=${CTRL_CN_DIR}"ctrl_denois_CN_NofCTRL="${N}"_PoN-"${CTRL_PON_SEX_SELECT}"_LOSS="${SEGM_MINLOSS}"_GAIN="${SEGM_MINGAIN}".tsv"
				CTRL_CN_FILE_SEGMENT_BASE=${CTRL_CN_DIR}"ctrl_denois_CN_NofCTRL="${N}"_PoN-"${CTRL_PON_SEX_SELECT}"_LOSS="${SEGM_MINLOSS}"_GAIN="${SEGM_MINGAIN}

				### use control denoised CN to prepare pre-file
				#============================================================================
				if ! ls ${CTRL_CN_FILE} 1> /dev/null 2>&1; then
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ preparing CN rds pre-file using N="${N}" controls - for PoN-"${CTRL_PON_SEX_SELECT}${COL2}
					echo -e ${COL1}${END}${COL2}
					mkdir -p ${CTRL_CN_DIR}
					Rscript --vanilla ${R}"ctrl_CN.R" \
						${INTERVALS_PROCESSED} \
						${CTRL_DENOIS_DIR} \
						${PON_SEX_PATTERN} \
						${N} \
						${CTRL_CN_FILE}
					if [[ ${CHECKPOINTS} == "yes" ]]; then 
						if ls "status_ok" 1> /dev/null 2>&1; then
							rm "status_ok"
						else
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: R script was halted"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Stage: CN in controls"
							echo -e "   File not generated: "${CTRL_CN_FILE}
							echo -e "   R script: ctrl_CN.R"
							echo -e ""${COL2}
							touch ./error
							exit
						fi
						if ! ls ${CTRL_CN_FILE} 1> /dev/null 2>&1; then
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Stage: CN in controls"
							echo -e "   File not generated: "${CTRL_CN_FILE}
							echo -e ""${COL2}
							touch ./error
							exit
						fi
					fi
				fi
				#============================================================================

				### process CN in controls using given segmentation conditions 
				#============================================================================
				#CHECK - CN CTRLS
				if ! ls ${CTRL_CN_FILE_SEGMENT} 1> /dev/null 2>&1; then
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ processing CN using N="${N}${COL2}
					echo -e ${COL1}"◼︎ segmenation conditions: LOSS "${SEGM_MINLOSS}" / GAIN "${SEGM_MINGAIN}${COL2}
					echo -e ${COL1}${END}${COL2}
					Rscript --vanilla ${R}"ctrl_CN_process_by_segmentation.R" \
						${CTRL_CN_FILE} \
						${SEGM_MINLOSS} \
						${SEGM_MINGAIN} \
						${CTRL_CN_FILE_SEGMENT}
					if [[ ${CHECKPOINTS} == "yes" ]]; then 
						if ls "status_ok" 1> /dev/null 2>&1; then
							rm "status_ok"
						else
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: R script was halted"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Stage: CN in controls"
							echo -e "   File not generated: "${CTRL_CN_FILE_SEGMENT}
							echo -e "   R script: ctrl_CN_process_by_segmentation.R"
							echo -e ""${COL2}
							touch ./error
							exit
						fi
						if ! ls ${CTRL_CN_FILE_SEGMENT} 1> /dev/null 2>&1; then
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Stage: CN in controls"
							echo -e "   File not generated: "${CTRL_CN_FILE_SEGMENT}
							echo -e ""${COL2}
							touch ./error
							exit
						fi
					fi

					#transform bedGraph to BigWig
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ processing CN using N="${N}${COL2}
					echo -e ${COL1}"◼︎ segmenation conditions: LOSS "${SEGM_MINLOSS}" / GAIN "${SEGM_MINGAIN}${COL2}
					echo -e ${COL1}"◼︎ transforming bedGraph to BigWig"${COL2}
					echo -e ${COL1}${END}${COL2}

					for file in ${CTRL_CN_FILE_SEGMENT_BASE}"_freqof"*; do
						${BEDGRAPHTOBIGWIG} ${file} ${CONTIGS_SIZES} ${file}".bw"
						rm ${file}
					done
				#CHECK - CN CTRLS
				fi
				#============================================================================

			#CHECK - CTRL CN REQUIRED
			fi

		#CHECK - CTRL YES	
		fi
	###====================================================================================

	### prepare gnomad SNP database (continuation) - if needed (exclude noisy SNP positions detected in controls)
	###====================================================================================

		#CHECK - CTRL YES 
		if [ ${CTRL} == "yes" ]; then

			### use control maf to prepare pre-file
			### two outputs: 1) denoised vcf -> collection allelic count, and 2) potentialy denoised rds -> report 
			#============================================================================
			if [ ${GNOMAD_SELECTION} == "capture_filter_and_ctrl_denois" ] || [ ${GNOMAD_SELECTION} == "af_filter_and_ctrl_denois" ]; then

				[ ${GNOMAD_SELECTION} == "capture_filter_and_ctrl_denois" ] && GNOMAD_SELECTED_IN=${GNOMAD_CAPTURE}
				[ ${GNOMAD_SELECTION} == "capture_filter_and_ctrl_denois" ] && GNOMAD_SELECTED_OUT=${GNOMAD_CAPTURE_DENOIS}
				[ ${GNOMAD_SELECTION} == "af_filter_and_ctrl_denois" ] && GNOMAD_SELECTED_IN=${GNOMAD_AF}
				[ ${GNOMAD_SELECTION} == "af_filter_and_ctrl_denois" ] && GNOMAD_SELECTED_OUT=${GNOMAD_AF_DENOIS}
				
				if ! ls ${GNOMAD_SELECTED_OUT} 1> /dev/null 2>&1; then
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ denoising gnomAD SNP database using N="${N}" controls"${COL2}
					echo -e ${COL1}"◼︎ excluding SNP position where more than "${AFPERC}" of controls have"${COL2}
					echo -e ${COL1}"◼︎ sequencing depth below "${AFDEPTH}" and/or"${COL2}
					echo -e ${COL1}"◼︎ not clear AA-0.0/AB-0.5/BB-1.0 genotype, with "${AFDIF}" deviation tolerance"${COL2}
					echo -e ${COL1}${END}${COL2}
					${BGZIP} -dc ${GNOMAD_SELECTED_IN} > ${GNOMAD_DIR}"temp1.vcf"
					#Find noisy SNPs
					Rscript --vanilla ${R}"ctrl_noisy_SNP.R" \
						${GNOMAD_DIR}"temp1.vcf" \
						${CTRL_MAF_DIR} \
						${GNOMAD_COUNT_PATTERN}".rds" \
						${GNOMAD_DIR}"temp2" \
						${AFDIF} \
						${AFDEPTH} \
						${AFPERC} \
						${N}
					if [[ ${CHECKPOINTS} == "yes" ]]; then
						if ls "status_ok" 1> /dev/null 2>&1; then
							rm "status_ok"
						else
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: R script was halted"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Stage: noisy SNPs in controls"
							echo -e "   File not generated: "${GNOMAD_DIR}"temp2"
							echo -e "   R script: ctrl_noisy_SNP.R"
							echo -e ""${COL2}
							touch ./error
							exit
						fi
						if ! ls ${GNOMAD_DIR}"temp2" 1> /dev/null 2>&1; then
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Stage: noisy SNPs in controls"
							echo -e "   File not generated: "${GNOMAD_DIR}"temp2"
							echo -e ""${COL2}
							touch ./error
							exit
						fi
					fi

					#Create new vcf with excluded noisy data
					grep "#" ${GNOMAD_DIR}"temp1.vcf" > ${GNOMAD_DIR}"temp1_head"
					cat ${GNOMAD_DIR}"temp1_head" ${GNOMAD_DIR}"temp2" > ${GNOMAD_DIR}"temp3.vcf"
					${BGZIP} -c ${GNOMAD_DIR}"temp3.vcf" > ${GNOMAD_SELECTED_OUT}
					${TABIX} -p vcf ${GNOMAD_SELECTED_OUT}
					rm ${GNOMAD_DIR}"temp"*

					if [[ ${CHECKPOINTS} == "yes" ]]; then
						if ! ls ${GNOMAD_SELECTED_OUT} 1> /dev/null 2>&1; then
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Stage: noisy SNPs in controls"
							echo -e "   File not generated: "${GNOMAD_SELECTED_OUT}
							echo -e ""${COL2}
							touch ./error
							exit
						fi
					fi
				fi
			fi
			#============================================================================

		#CHECK - CTRL YES 
		fi

	# remove temp files
	SEGMENTATION_TEMP=${SEGMENTATION}"_temp.txt"
	if ls ${SEGMENTATION_TEMP} 1> /dev/null 2>&1; then
		rm ${SEGMENTATION_TEMP}
	fi
	###====================================================================================

#LOOP - LINE oF PROJECTS
done
#============================================================================

if [[ ${CHECKPOINTS} == "yes" ]]; then 
	if ls ./error 1> /dev/null 2>&1; then
		echo -e ${COL1}${BEG}${COL2}
		echo -e ${COL1}"◼︎ $(date)"${COL2}
		echo -e ${COL1_WARN}"◼︎ CHECKPOINT ERROR - error(s) detected during preparation and/or control processing (listed above)"
		echo -e "◼︎ Please resolve before processing."
		echo -e ${COL1}${END}${COL2}
		exit
	fi
fi

#just when CHECKPOINTS are off, this file is still generated
if ls "status_ok" 1> /dev/null 2>&1; then rm "status_ok"; fi

echo -e ${COL1}${BEG}${COL2}
echo -e ${COL1}"◼︎ $(date)"${COL2}
echo -e ${COL1}"◼︎ controls preparation finished"${COL2}
echo -e ${COL1}${END}${COL2}

####============================================================================
