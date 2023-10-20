#!/bin/bash

#prepare reference bed file
if ! ls ${INTERVALS_PROCESSED}"_target.bed" 1> /dev/null 2>&1; then

	echo -e ${COL1}${BEG}${COL2}
	echo -e ${COL1}"◼︎ $(date)"${COL2}
	echo -e ${COL1}"◼︎ PREPROCESSING - CNVkit: reference target preparation"${COL2}
	echo -e ${COL1}${END}${COL2}

	Rscript --vanilla ${R}"cnvkit_ref_target.R" \
		${INTERVALS_PROCESSED} \
		${INTERVALS_PROCESSED}"_target.bed"

	#create empty antitarget
	touch ${INTERVALS_PROCESSED}"_antitarget.bed"

	if [[ ${CHECKPOINTS} == "yes" ]]; then 
		if ls "status_ok" 1> /dev/null 2>&1; then
			rm "status_ok"
		else
			echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: R script was halted"
			echo -e "   Project: "${PROJECT_ID}
			echo -e "   Stage: CNVkit - target preparation"
			echo -e "   File not generated: "${INTERVALS_PROCESSED}"_target.bed"
			echo -e "   R script: cnvkit_ref_target.R"
			echo -e ""${COL2}
			touch ./error
			exit
		fi
		if ! ls ${INTERVALS_PROCESSED}"_target.bed" 1> /dev/null 2>&1; then
			echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
			echo -e "   Project: "${PROJECT_ID}
			echo -e "   Stage: CNVkit - target preparation"
			echo -e "   File not generated: "${INTERVALS_PROCESSED}"_target.bed"
			echo -e ""${COL2}
			touch ./error
			exit
		fi
	fi
fi



#Collect coverage for each control
###====================================================================================

#CHECK - CTRL YES 
if [[ ${CTRL} == "yes" ]]; then

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

		[[ ${WAY_TO_BAM} == "find_in_dir" ]] && BAMFILE=`find ${CTRL_BAM_DIR} -type f -name ${CTRL_ID}*${CTRL_BAM_PATTERN}`
		[[ ${WAY_TO_BAM} == "absolute" ]] && BAMFILE=${CTRL_BAM_DIR}${CTRL_BAM_PATH}


#-----CNVkit active code

		if [[ ${CTRL_SEX} == "unk" ]]; then 
			if ls ${CTRL_COUNT_CNN_DIR}${MAIN_ID}"_"${CTRL_ID}"_M_targetcoverage.cnn" 1> /dev/null 2>&1; then CTRL_SEX="M"; fi
			if ls ${CTRL_COUNT_CNN_DIR}${MAIN_ID}"_"${CTRL_ID}"_F_targetcoverage.cnn" 1> /dev/null 2>&1; then CTRL_SEX="F"; fi
			if ls ${CTRL_COUNT_CNN_DIR}${MAIN_ID}"_"${CTRL_ID}"_unk_targetcoverage.cnn" 1> /dev/null 2>&1; then CTRL_SEX="unk"; fi
		fi

		RESULTS_DATA_PATTERN=${MAIN_ID}"_"${CTRL_ID}"_"${CTRL_SEX}
		ECHO_PATTERN=${MAIN_ID}" - "${CTRL_ID}



		CTRL_COUNT_DIR=${CTRL_COUNT_CNN_DIR}
		RESULTS_DATA_PATTERN=${MAIN_ID}"_"${CTRL_ID}"_"${CTRL_SEX}
		ECHO_PATTERN=${MAIN_ID}" - "${CTRL_ID}

		OUTPUT=${CTRL_COUNT_DIR}${RESULTS_DATA_PATTERN}"_targetcoverage.cnn"
		OUTPUT_ANTI=${CTRL_COUNT_DIR}${RESULTS_DATA_PATTERN}"_antitargetcoverage.cnn"
		if ! ls ${OUTPUT} 1> /dev/null 2>&1; then

			echo -e ${COL1}${BEG}${COL2}
			echo -e ${COL1}"◼︎ $(date)"${COL2}
			echo -e ${COL1}"◼︎ controls prep - "${ECHO_PATTERN}": CNVkit collecting coverage"${COL2}
			echo -e ${COL1}${END}${COL2}
			mkdir -p ${CTRL_COUNT_DIR}

			source ${CONDA}
			conda activate ${CNVKIT_ENV}

			cnvkit.py coverage ${BAMFILE} ${INTERVALS_PROCESSED}"_target.bed" -o ${OUTPUT}
			cnvkit.py coverage ${BAMFILE} ${INTERVALS_PROCESSED}"_antitarget.bed" -o ${OUTPUT_ANTI}

			conda deactivate

			if [[ ${CHECKPOINTS} == "yes" ]]; then 
				if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
					echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
					echo -e "   Project: "${PROJECT_ID}
					echo -e "   Main ID: "${MAIN_ID}
					echo -e "   Control ID: "${CTRL_ID}
					echo -e "   Stage: CNVkit - collection coverage"
					echo -e "   File not generated: "${OUTPUT}
					echo -e ""${COL2}
					touch ./error
					exit
				fi
			fi




###WORKING PART - CHANGE TO CNN prediction of sex

			#CHECK SEX
			if [[ ${CTRL_SEX_TEST} == "default" || ${CTRL_SEX_TEST} == "yes" ]]; then
				echo -e ${COL1}${BEG}${COL2}
				echo -e ${COL1}"◼︎ $(date)"${COL2}
				echo -e ${COL1}"◼︎ controls prep - "${ECHO_PATTERN}": CNVkit - predicting sex"${COL2}
				echo -e ${COL1}${END}${COL2}
				CTRL_SEX_PREDICTION=`Rscript --vanilla ${R}"ctrl_predict_sex.R" ${CTRL_COUNT_DIR}${RESULTS_DATA_PATTERN}"_targetcoverage.cnn"`
				
				if [[ ${CHECKPOINTS} == "yes" ]]; then 
					if ls "status_ok" 1> /dev/null 2>&1; then
						rm "status_ok"
					else
						echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: R script was halted"
						echo -e "   Project: "${PROJECT_ID}
						echo -e "   Main ID: "${MAIN_ID}
						echo -e "   Control ID: "${CTRL_ID}
						echo -e "   Stage: CNVkit - sex prediction"
						echo -e "   R script: ctrl_predict_sex.R"
						echo -e ""${COL2}
						touch ./error
						exit
					fi

					if [[ ${CTRL_SEX_PREDICTION} != "M" && ${CTRL_SEX_PREDICTION} != "F" && ${CTRL_SEX_PREDICTION} != "unk" ]]; then
						echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: sex  predicted with unexpected value"
						echo -e "   Project: "${PROJECT_ID}
						echo -e "   Main ID: "${MAIN_ID}
						echo -e "   Control ID: "${CTRL_ID}
						echo -e "   Stage: CNVkit - sex prediction"
						echo -e "   Sex expected to be F, M or unk but predicted as: "${CTRL_SEX_PREDICTION}
						echo -e ""${COL2}
						touch ./error
						exit
					fi
				fi

				#change name of the file by the predicted
				if [[ ${CTRL_SEX} == "unk" ]]; then 
					if [[ ${CTRL_SEX_PREDICTION} != "unk" ]]; then
						echo -e ${COL1_WARN}"◼︎ INFO CNVkit - Sex of control with unknonw sex predicted"
						echo -e "   Project: "${PROJECT_ID}
						echo -e "   Main ID: "${MAIN_ID}
						echo -e "   Control ID: "${CTRL_ID}
						echo -e "   Given sex by user: "${CTRL_SEX}
						echo -e "   Predicted sex: "${CTRL_SEX_PREDICTION}
						echo -e ""${COL2}
						mv ${CTRL_COUNT_DIR}${RESULTS_DATA_PATTERN}"_targetcoverage.cnn" ${CTRL_COUNT_DIR}${MAIN_ID}"_"${CTRL_ID}"_"${CTRL_SEX_PREDICTION}"_targetcoverage.cnn"
						mv ${CTRL_COUNT_DIR}${RESULTS_DATA_PATTERN}"_antitargetcoverage.cnn" ${CTRL_COUNT_DIR}${MAIN_ID}"_"${CTRL_ID}"_"${CTRL_SEX_PREDICTION}"_antitargetcoverage.cnn"
					fi
				fi

				#if sex is defined in the table, check if it looks correct and exit if not (as long as it was possible to predict)
				if [[ ${CTRL_SEX} != "unk" ]]; then 
					if [[ ${CTRL_SEX_PREDICTION} != "unk" ]]; then 
						if [[ ${CTRL_SEX} == ${CTRL_SEX_PREDICTION} ]]; then 
							echo -e ${COL1_WARN}"◼︎ INFO CNVkit - Control sex successfully verified"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Main ID: "${MAIN_ID}
							echo -e "   Control ID: "${CTRL_ID}
							echo -e "   Sex given in the table: "${CTRL_SEX}
							echo -e "   Sex predicted by coverage: "${CTRL_SEX_PREDICTION}
							echo -e ""${COL2}
						fi
						if [[ ${CTRL_SEX} != ${CTRL_SEX_PREDICTION} ]]; then 
							echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: !!! CHECK SEX OF CONTROLS !!!"
							echo -e "   Project: "${PROJECT_ID}
							echo -e "   Main ID: "${MAIN_ID}
							echo -e "   Control ID: "${CTRL_ID}
							echo -e "   Stage: CNVkit - sex test"
							echo -e "   Sex given in the controls master is not matching sex predicted by coverage"
							echo -e "   Sex given in the table: "${CTRL_SEX}
							echo -e "   Sex predicted by coverage: "${CTRL_SEX_PREDICTION}
							echo -e ""${COL2}
							touch ./error
							exit
						fi
					fi
				fi

			#CHECK SEX
			fi


###WORKING PART - CHANGE TO CNN prediction of sex








		fi


		# #double chrX counts in males
		# if [[ ${AUTOSOME_MODE} != "A" && ${CTRL_SEX} == "M" ]]; then
		# 	OUTPUT=${CTRL_COUNT_DIR}${RESULTS_DATA_PATTERN}"_targetcoverage_doublegon.cnn"
		# 	if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
		# 		echo -e ${COL1}${BEG}${COL2}
		# 		echo -e ${COL1}"◼︎ $(date)"${COL2}
		# 		echo -e ${COL1}"◼︎ controls prep - "${ECHO_PATTERN}": CNVkit - doubling gonosomal coverage for this male control"${COL2}
		# 		echo -e ${COL1}${END}${COL2}
		# 		Rscript --vanilla ${R}"cnvkit_gonosome_double.R" \
		# 			${CTRL_COUNT_DIR}${RESULTS_DATA_PATTERN}"_targetcoverage.cnn" \
		# 			${OUTPUT}

		# 		if [[ ${CHECKPOINTS} == "yes" ]]; then 
		# 			if ls "status_ok" 1> /dev/null 2>&1; then
		# 				rm "status_ok"
		# 			else
		# 				echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: R script was halted"
		# 				echo -e "   Project: "${PROJECT_ID}
		# 				echo -e "   Main ID: "${MAIN_ID}
		# 				echo -e "   Control ID: "${CTRL_ID}
		# 				echo -e "   Stage: CNVkit - doubling gonosomal counts for male control"
		# 				echo -e "   File not generated: "${OUTPUT}
		# 				echo -e "   R script: cnvkit_gonosome_double.R"
		# 				echo -e ""${COL2}
		# 				touch ./error
		# 				exit
		# 			fi
		# 			if ! ls ${OUTPUT} 1> /dev/null 2>&1; then
		# 				echo -e ${COL1_INFO}"◼︎ CHECKPOINT ERROR - CONTROL PREP: file not generated"
		# 				echo -e "   Project: "${PROJECT_ID}
		# 				echo -e "   Main ID: "${MAIN_ID}
		# 				echo -e "   Control ID: "${CTRL_ID}
		# 				echo -e "   Stage: CNVkit - doubling gonosomal counts for male control"
		# 				echo -e "   File not generated: "${OUTPUT}
		# 				echo -e ""${COL2}
		# 				touch ./error
		# 				exit
		# 			fi
		# 		fi
		# 	fi
		# #double gonosomal counts in males
		# fi

#-----CNVkit active code





	#LOOP - CTRLS
	done

#CHECK - CTRL YES 
fi
###====================================================================================








