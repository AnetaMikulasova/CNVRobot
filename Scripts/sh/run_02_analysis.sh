#!/bin/bash

source ./setting.sh

echo -e ${COL1}${BEG}${COL2}
echo -e ${COL1}"◼︎ $(date)"${COL2}
echo -e ${COL1}"◼︎ samples analysis ..."${COL2}
echo -e ${COL1}${END}${COL2}

#processing samples
#====================================================================================

#LOOP - LINE oF PROJECTS
#============================================================================
cat ${MASTERPROJECTS} | awk -F"\t" '$1=="yes"' | while read -r line || [[ -n "$line" ]]; do

	PROJECT_ID=`extract_from_master 'PROJECT_ID' ${MASTERPROJECTS}`
	PROJECT_TYPE=`extract_from_master 'PROJECT_TYPE' ${MASTERPROJECTS}`
	SEQ_TYPE=`extract_from_master 'SEQ_TYPE' ${MASTERPROJECTS}`
	CAPTURE_ID=`extract_from_master 'CAPTURE_ID' ${MASTERPROJECTS}`
	CAPTURE_FILE=`extract_from_master 'CAPTURE_FILE' ${MASTERPROJECTS}`
	GENOME_VERSION=`extract_from_master 'GENOME_VERSION' ${MASTERPROJECTS}`
	REF=`extract_from_master 'REF' ${MASTERPROJECTS}`
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


	#count number of controls
	Nx=`cat ${MASTERCONTROLS} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION)print}' | wc -l`
	# Nx=`ls ${CTRL_DENOIS_DIR}* | grep "denoisedCR_PoN" | wc -l`
	[ ${CTRL} == "yes" ] && N=$(echo $Nx | tr -d ' ')
	[ ${CTRL} == "no" ] && N=0

	SUPPORT_FILES=${SUPPORT_FILES_DIR}"/"${PROJECT_ID}_${CAPTURE_ID}"_"${GENOME_VERSION}"_bin"${BIN}"bp_pad"${PADDING}"bp/"
	CTRL_COUNT_COUNTS_DIR=${SUPPORT_FILES}"Controls/counts/"
	CTRL_COUNT_HDF_DIR=${SUPPORT_FILES}"Controls/counts/hdf5/"
    CTRL_COUNT_MAF_DIR=${SUPPORT_FILES}"Controls/maf/"

	SAMPLE_PoN_DIR=${SUPPORT_FILES}"Samples/PoN/NofCTRL="${N}"/"${PROJECT_ID}"/"
	if [ ${CTRL} == "yes" ]; then 
		mkdir -p ${SAMPLE_PoN_DIR}
	fi

	INTERVALS=${SUPPORT_FILES}"Capture_Intervals/"
	INTERVALS_PROCESSED=${INTERVALS}${CAPTURE_ID}"_preprocessed_intervals.interval_list"
	INTERVALS_ANNOTATED=${INTERVALS}${CAPTURE_ID}"_annotated_intervals.interval_list"
	CONTIGS_SIZES=${INTERVALS}${GENOME_VERSION}"_contigs_sizes.txt"


	#define segmentation conditions
	###====================================================================================
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
	#recognize if user gives segmentation
	fi

	### calculate segmentation conditions when smart option
	# requires annotated intervals
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
	fi
	###====================================================================================


	#assign databases automatically
	###====================================================================================
	if [ ${GENOME_VERSION} == "GRCh37-hg19" ]; then
		if [ `grep "chr" ${CONTIGS_SIZES} | wc -l` == 0 ]; then
			PLOT_ASSEMBLY="hg19"
			GNOMAD=${DATABASES_DIR}"GNOMAD/gnomad_2.1_GRCh37/gnomad_2.1_combine_processed_GRCh37_22XY.vcf.gz"
			DGV=${DATABASES_DIR}"DGV/DGV_GRCh37_processed_variants_2020-02-25.rds"
			UCSC_CHROMOSOME_BANDS=${DATABASES_DIR}"UCSC/chr_bands/GRCh37_cytoBand.txt"
			UCSC_CHROMOSOME_CENTROMERE=${DATABASES_DIR}"UCSC/chr_centromere/GRCh37_centromeres.txt"
			MAPPABILITY=${DATABASES_DIR}"MAPPABILITY/mappability_GRCh37.rds"
		else
			PLOT_ASSEMBLY="hg19"
			GNOMAD=${DATABASES_DIR}"GNOMAD/gnomad_2.1_hg19_liftover/gnomad_2.1_combine_processed_hg19_22XY.vcf.gz"
			DGV=${DATABASES_DIR}"DGV/DGV_hg19_processed_variants_2020-02-25.rds"
			UCSC_CHROMOSOME_BANDS=${DATABASES_DIR}"UCSC/chr_bands/hg19_cytoBand.txt"
			UCSC_CHROMOSOME_CENTROMERE=${DATABASES_DIR}"UCSC/chr_centromere/hg19_centromeres.txt"
			MAPPABILITY=${DATABASES_DIR}"MAPPABILITY/mappability_hg19.rds"
		fi
		if [ ${GENE_ANNOTATION} == "default" ]; then
			GENE_ANNOTATION=${DATABASES_DIR}"/GENE_ANNOTATION/RefSeqAll_hg19_191206.txt"
		fi
	fi
	if [ ${GENOME_VERSION} == "GRCh38-hg38" ]; then
		if [ `grep "chr" ${CONTIGS_SIZES} | wc -l` != 0 ]; then
			PLOT_ASSEMBLY="hg38"
	 		GNOMAD=${DATABASES_DIR}"GNOMAD/gnomad_2.1_hg38_liftover/gnomad_2.1_combine_processed_hg38_22XY.vcf.gz"
	 		DGV=${DATABASES_DIR}"DGV/DGV_GRCh38hg38_processed_variants_2020-02-25.rds"
	 		UCSC_CHROMOSOME_BANDS=${DATABASES_DIR}"UCSC/chr_bands/GRCh38_hg38_cytoBand.txt"
	 		UCSC_CHROMOSOME_CENTROMERE=${DATABASES_DIR}"UCSC/chr_centromere/GRCh38_hg38_centromeres.txt"
	 		MAPPABILITY=${DATABASES_DIR}"MAPPABILITY/mappability_hg38GRCh38.rds"
	 	else
			PLOT_ASSEMBLY="hg38"
	 		GNOMAD=${DATABASES_DIR}"GNOMAD/gnomad_2.1_GRCh38_liftover/gnomad_2.1_combine_processed_GRCh38_22XY.vcf.gz"
	 		DGV=${DATABASES_DIR}"DGV/DGV_GRCh38_processed_variants_2020-02-25.rds"
	 		UCSC_CHROMOSOME_BANDS=${DATABASES_DIR}"UCSC/chr_bands/GRCh38_cytoBand.txt"
	 		UCSC_CHROMOSOME_CENTROMERE=${DATABASES_DIR}"UCSC/chr_centromere/GRCh38_centromeres.txt"
	 		MAPPABILITY=${DATABASES_DIR}"MAPPABILITY/mappability_GRCh38.rds"
	 	fi
	 	if [ ${GENE_ANNOTATION} == "default" ]; then
			GENE_ANNOTATION=${DATABASES_DIR}"/GENE_ANNOTATION/RefSeqAll_hg38_191206.txt"
		fi
	fi
	###====================================================================================


	#process genes database if required processed file is not found
	#============================================================================
	#CHECK - GENES PROCESSED
	GENE_ANNOTATION_PROCESSED="${GENE_ANNOTATION%.*}"
	GENE_ANNOTATION_PROCESSED=${GENE_ANNOTATION_PROCESSED}"_processed.txt"
	if ! ls ${GENE_ANNOTATION_PROCESSED} 1> /dev/null 2>&1; then
		echo -e ${COL1}${BEG}${COL2}
		echo -e ${COL1}"◼︎ $(date)"${COL2}
		echo -e ${COL1}"◼︎ PREPROCESSING - processed gene database not found, processing original source (needed just once)"${COL2}
		echo -e ${COL1}${END}${COL2}
		#process gene database
		Rscript --vanilla ${R}"process_gene_database.R" \
			${GENE_ANNOTATION} \
			${GENE_ANNOTATION_PROCESSED}
		#CHECK - GENES PROCESSED
	fi
	#============================================================================


	### assign variables
	#============================================================================
	GNOMAD_ID=$(basename ${GNOMAD})
	GNOMAD_ID=`echo "${GNOMAD_ID//.vcf/}"`
	GNOMAD_ID=`echo "${GNOMAD_ID//.gz/}"`
	# GNOMAD_ID="${GNOMAD_ID%.*}"
	# GNOMAD_ID="${GNOMAD_ID%.*}"
	GNOMAD_DIR=${SUPPORT_FILES}"gnomAD/"

	GNOMAD_CAPTURE=${GNOMAD_DIR}${GNOMAD_ID}"_1_filtered_for_capture.vcf.gz"
	GNOMAD_CAPTURE_DENOIS=${GNOMAD_DIR}${GNOMAD_ID}"_2_filter_for_capture_and_denoised_NofCTRL="${N}".vcf.gz"
	# AF_GNOMAD="0.001"
	GNOMAD_AF=${GNOMAD_DIR}${GNOMAD_ID}"_1_filtered_for_af_"${AF_GNOMAD}".vcf.gz"
	GNOMAD_AF_DENOIS=${GNOMAD_DIR}${GNOMAD_ID}"_2_filtered_for_af_"${AF_GNOMAD}"_and_denoised_NofCTRL="${N}".vcf.gz"

	[ ${GNOMAD_SELECTION} == "capture_filter" ] && GNOMAD_SELECTED=${GNOMAD_CAPTURE} && GNOMAD_COUNT_PATTERN="_filter_capture"
	[ ${GNOMAD_SELECTION} == "capture_filter_and_ctrl_denois" ]  && GNOMAD_SELECTED=${GNOMAD_CAPTURE_DENOIS} && GNOMAD_COUNT_PATTERN="_filter_capture_denois"
	[ ${GNOMAD_SELECTION} == "af_filter" ] && GNOMAD_SELECTED=${GNOMAD_AF} && GNOMAD_COUNT_PATTERN="_filter_af_"${AF_GNOMAD}
	[ ${GNOMAD_SELECTION} == "af_filter_and_ctrl_denois" ] && GNOMAD_SELECTED=${GNOMAD_AF_DENOIS} && GNOMAD_COUNT_PATTERN="_filter_af_"${AF_GNOMAD}"_denois"

	CTRL_CN_DIR=${SUPPORT_FILES}"Controls/CN/"
	CTRL_CN_FILE=${CTRL_CN_DIR}"ctrl_denois_CN_NofCTRL="${N}"_PoN-"${CTRL_PON_SEX_SELECT}".rds"
	CTRL_CN_FILE_SEGMENT=${CTRL_CN_DIR}"ctrl_denois_CN_NofCTRL="${N}"_PoN-"${CTRL_PON_SEX_SELECT}"_LOSS="${SEGM_MINLOSS}"_GAIN="${SEGM_MINGAIN}".tsv"

	#============================================================================

	#perform test if number of processed controls is the same as expected
	#====================================================================================
	[ ${CTRL} == "yes" ] && Ntestx=`ls ${CTRL_COUNT_HDF_DIR}* | grep "_counts.hdf5" | wc -l`
	[ ${CTRL} == "no" ] && Ntestx=0
	Ntest=$(echo ${Ntestx} | tr -d ' ')
	if [ ${N} != ${Ntest} ]; then
		echo -e ${COL1}${BEG}${COL2}
		echo -e ${COL1}"◼︎ $(date)"${COL2}
		echo -e ${COL1_WARN}"◼︎ !!! CHECK CONTROLS, NOT EXPECTED NUMBER OF THEM !!!"${COL2}
		echo -e ${COL1}${END}${COL2}
		exit
	fi
	#====================================================================================


	#LOOP - LINE oF MASTER
	#============================================================================
	cat ${MASTERSAMPLES} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION)print}' | while read -r line || [[ -n "$line" ]]; do

	 	MAIN_ID=`extract_from_master 'MAIN_ID' ${MASTERSAMPLES}`
	 	SAMPLE_1_ID=`extract_from_master 'SAMPLE1_ID' ${MASTERSAMPLES}`
	 	SAMPLE_1_TYPE=`extract_from_master 'SAMPLE1_TYPE' ${MASTERSAMPLES}`
	 	SAMPLE_1_SEX=`extract_from_master 'SAMPLE1_SEX' ${MASTERSAMPLES}`
	 	SAMPLE_1_BAM_PATH=`extract_from_master 'SAMPLE1_PATH_TO_BAM' ${MASTERSAMPLES}`
	 	SAMPLE_2_ID=`extract_from_master 'SAMPLE2_ID' ${MASTERSAMPLES}`
	 	SAMPLE_2_TYPE=`extract_from_master 'SAMPLE2_TYPE' ${MASTERSAMPLES}`
	 	SAMPLE_2_SEX=`extract_from_master 'SAMPLE2_SEX' ${MASTERSAMPLES}`
	 	SAMPLE_2_BAM_PATH=`extract_from_master 'SAMPLE2_PATH_TO_BAM' ${MASTERSAMPLES}`
	 	SAMPLE_3_ID=`extract_from_master 'SAMPLE3_ID' ${MASTERSAMPLES}`
	 	SAMPLE_3_TYPE=`extract_from_master 'SAMPLE3_TYPE' ${MASTERSAMPLES}`
	 	SAMPLE_3_SEX=`extract_from_master 'SAMPLE3_SEX' ${MASTERSAMPLES}`
	 	SAMPLE_3_BAM_PATH=`extract_from_master 'SAMPLE3_PATH_TO_BAM' ${MASTERSAMPLES}`

	 	#change variables that are empty to na
	 	[[ ${SAMPLE_1_ID} == "" ]] && SAMPLE_1_ID="NA"
	 	[[ ${SAMPLE_1_TYPE} == "" ]] && SAMPLE_1_TYPE="NA"
	 	[[ ${SAMPLE_1_SEX} == "" ]] && SAMPLE_1_SEX="NA"
	 	[[ ${SAMPLE_1_BAM_PATH} == "" ]] && SAMPLE_1_BAM_PATH="NA"
	 	[[ ${SAMPLE_2_ID} == "" ]] && SAMPLE_2_ID="NA"
	 	[[ ${SAMPLE_2_TYPE} == "" ]] && SAMPLE_2_TYPE="NA"
	 	[[ ${SAMPLE_2_SEX} == "" ]] && SAMPLE_2_SEX="NA"
	 	[[ ${SAMPLE_2_BAM_PATH} == "" ]] && SAMPLE_2_BAM_PATH="NA"
	 	[[ ${SAMPLE_3_ID} == "" ]] && SAMPLE_3_ID="NA"
	 	[[ ${SAMPLE_3_TYPE} == "" ]] && SAMPLE_3_TYPE="NA"
	 	[[ ${SAMPLE_3_SEX} == "" ]] && SAMPLE_3_SEX="NA"
	 	[[ ${SAMPLE_3_BAM_PATH} == "" ]] && SAMPLE_3_BAM_PATH="NA"

		RESULTS_DIR_GO=${RESULTS_DIR}${PROJECT_ID}"/"${PROJECT_ID}"_"${CAPTURE_ID}"_"${GENOME_VERSION}"_bin"${BIN}"bp_pad"${PADDING}"bp_NofCTRLS="${N}"/"${MAIN_ID}"/"
		RESULTS_DATA_DIR=${RESULTS_DIR_GO}"data/"
		mkdir -p ${RESULTS_DATA_DIR}
		RESULTS_PLOT_DIR=${RESULTS_DIR_GO}"plot_"${SEGMENTATIONCONDITIONS}"/"
		RESULTS_IGV_BED_DIR=${RESULTS_DIR_GO}"IGV_"${SEGMENTATIONCONDITIONS}"/"

		#process samples: read counts, prepare PoN
		#====================================================================================
		#LOOP - SAMPLES PROCESSING
		for CASE_ID in ${SAMPLE_1_ID} ${SAMPLE_2_ID} ${SAMPLE_3_ID}; do

			[ ${CASE_ID} == ${SAMPLE_1_ID} ] && CASE_TYPE=${SAMPLE_1_TYPE}
			[ ${CASE_ID} == ${SAMPLE_2_ID} ] && CASE_TYPE=${SAMPLE_2_TYPE}
			[ ${CASE_ID} == ${SAMPLE_3_ID} ] && CASE_TYPE=${SAMPLE_3_TYPE}

			[ ${CASE_ID} == ${SAMPLE_1_ID} ] && CASE_SEX=${SAMPLE_1_SEX}
			[ ${CASE_ID} == ${SAMPLE_2_ID} ] && CASE_SEX=${SAMPLE_2_SEX}
			[ ${CASE_ID} == ${SAMPLE_3_ID} ] && CASE_SEX=${SAMPLE_3_SEX}

			[ ${CASE_ID} == ${SAMPLE_1_ID} ] && CASE_BAM_PATH=${SAMPLE_1_BAM_PATH}
			[ ${CASE_ID} == ${SAMPLE_2_ID} ] && CASE_BAM_PATH=${SAMPLE_2_BAM_PATH}
			[ ${CASE_ID} == ${SAMPLE_3_ID} ] && CASE_BAM_PATH=${SAMPLE_3_BAM_PATH}

			[ ${WAY_TO_BAM} == "find_in_dir" ] && BAMFILE=`find ${SMPL_BAM_DIR} -type f -name ${CASE_ID}*${SMPL_BAM_PATTERN}`
			[ ${WAY_TO_BAM} == "absolute" ] && BAMFILE=${SMPL_BAM_DIR}${CASE_BAM_PATH}

			#CHECK - SAMPLES EXISTS
			if ls ${BAMFILE} 1> /dev/null 2>&1; then

				if [ ${WAY_TO_BAM} == "find_in_dir" ]; then
					#test n of bam files, has to be one! - if "find_in_dir"
					BAMFILE_N=$(echo ${BAMFILE} | wc -w | tr -d ' ')
					if [ $BAMFILE_N != 1 ]; then
						echo -e ${COL1}"\n""◼︎ ERROR - unexpected # of bam file(s) of control:"${COL2}
						echo -e ${COL1}"◼︎ ERROR - "${CASE_ID}" - N="${BAMFILE_N}"\n"${COL2}
						exit
					fi
				fi

				RESULTS_DATA_PATTERN=${RESULTS_DATA_DIR}${CASE_TYPE}"_"${CASE_ID}"_"${CASE_SEX}

				#colect read counts
				#============================================================================
				#LOOP - READ COUNTS FORMAT
				for FILE_TYPE in "HDF5" "TSV"; do
					[ ${FILE_TYPE} == "HDF5" ] && FILE_TYPE_END="hdf5"
					[ ${FILE_TYPE} == "TSV" ] && FILE_TYPE_END="tsv"
					#CHECK - READ COUNTS
					if ! ls ${RESULTS_DATA_PATTERN}"_counts."${FILE_TYPE_END} 1> /dev/null 2>&1; then

						#CHECK - AS CTRL EXISTS
						if ls ${CTRL_COUNT_COUNTS_DIR}${FILE_TYPE_END}"/"${MAIN_ID}"_"${CASE_ID}"_"${CASE_SEX}"_counts."${FILE_TYPE_END} 1> /dev/null 2>&1; then
							echo -e ${COL1}${BEG}${COL2}
							echo -e ${COL1}"◼︎ $(date)"${COL2}
							echo -e ${COL1}"◼︎ "${MAIN_ID}" - "${CASE_TYPE}" - "${CASE_ID}": collecting read count - "${FILE_TYPE_END}${COL2} 
							echo -e ${COL1}"◼︎ collected as control, copying data..."${COL2} 							
							echo -e ${COL1}${END}${COL2}
							cp ${CTRL_COUNT_COUNTS_DIR}${FILE_TYPE_END}"/"${MAIN_ID}"_"${CASE_ID}"_"${CASE_SEX}"_counts."${FILE_TYPE_END} ${RESULTS_DATA_PATTERN}"_counts."${FILE_TYPE_END}
						else
						echo -e ${COL1}${BEG}${COL2}
						echo -e ${COL1}"◼︎ $(date)"${COL2}
						echo -e ${COL1}"◼︎ "${MAIN_ID}" - "${CASE_TYPE}" - "${CASE_ID}": collecting read count - "${FILE_TYPE_END}${COL2} 
						echo -e ${COL1}${END}${COL2}
						${GATK} CollectReadCounts \
						    -R ${REF} \
						    -I ${BAMFILE} \
						    -L ${INTERVALS_PROCESSED} \
						    --interval-merging-rule OVERLAPPING_ONLY \
						    --format ${FILE_TYPE} \
						    -O ${RESULTS_DATA_PATTERN}"_counts."${FILE_TYPE_END}
						#CHECK - AS CTRL EXISTS
						fi
					#CHECK - READ COUNTS
					fi
				#LOOP - READ COUNTS FORMAT
				done
				#============================================================================

				#collect allelic counts and convert to rds (-RF GoodCigarReadFilter \)
				#============================================================================
				#CHECK - MAF
				if ! ls ${RESULTS_DATA_PATTERN}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".rds" 1> /dev/null 2>&1; then

					#CHECK - AS CTRL EXISTS
					if ls ${CTRL_COUNT_MAF_DIR}${MAIN_ID}"_"${CASE_ID}"_"${CASE_SEX}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".rds" 1> /dev/null 2>&1; then
						echo -e ${COL1}${BEG}${COL2}
						echo -e ${COL1}"◼︎ $(date)"${COL2}
						echo -e ${COL1}"◼︎ "${MAIN_ID}" - "${CASE_TYPE}" - "${CASE_ID}": collection allelic counts"${COL2}
						echo -e ${COL1}"◼︎ collected as control, copying data..."${COL2} 	
						echo -e ${COL1}${END}${COL2}
						cp ${CTRL_COUNT_MAF_DIR}${MAIN_ID}"_"${CASE_ID}"_"${CASE_SEX}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".rds" ${RESULTS_DATA_PATTERN}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".rds"
					else
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ "${MAIN_ID}" - "${CASE_TYPE}" - "${CASE_ID}": collection allelic counts"${COL2} 
					echo -e ${COL1}${END}${COL2}

					${GATK} CollectAllelicCounts \
					   	-L ${GNOMAD_SELECTED} \
					   	-I ${BAMFILE} \
					   	-R ${REF} \
					   	-O ${RESULTS_DATA_PATTERN}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".tsv"
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ "${MAIN_ID}" - "${CASE_TYPE}" - "${CASE_ID}": converting allelic counts file to rds"${COL2} 
					echo -e ${COL1}${END}${COL2}	
					Rscript --vanilla ${R}convert_alleliccounts_to_rds.R \
							${RESULTS_DATA_PATTERN}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".tsv" \
							${RESULTS_DATA_PATTERN}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".rds"
					rm ${RESULTS_DATA_PATTERN}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".tsv"
					#CHECK - AS CTRL EXISTS
					fi
				#CHECK - MAF
				fi
				#============================================================================

				#run futher processing for all required PoN
				#============================================================================

				[ ${SMPL_PON_SEX_SELECT} == "M" ] && PON_SEX1="M" && PON_SEX2="NA" && PON_SEX3="NA" && PON_SEX4="NA " && PON_SEX5="NA"
				[ ${SMPL_PON_SEX_SELECT} == "F" ] && PON_SEX1="F" && PON_SEX2="NA" && PON_SEX3="NA" && PON_SEX4="NA" && PON_SEX5="NA"
				[ ${SMPL_PON_SEX_SELECT} == "mixed" ] && PON_SEX1="mixed" && PON_SEX2="NA" && PON_SEX3="NA" && PON_SEX4="NA" && PON_SEX5="NA"
				[ ${SMPL_PON_SEX_SELECT} == "both_separated" ] && PON_SEX1="M" && PON_SEX2="F" && PON_SEX3="NA" && PON_SEX4="NA" && PON_SEX5="NA"
				[ ${SMPL_PON_SEX_SELECT} == "matched_main" ] && PON_SEX1="matched_main" && PON_SEX2="NA" && PON_SEX3="NA" && PON_SEX4="NA" && PON_SEX5="NA"
				[ ${SMPL_PON_SEX_SELECT} == "matched_each" ] && PON_SEX1="matched_each" && PON_SEX2="NA" && PON_SEX3="NA" && PON_SEX4="NA" && PON_SEX5="NA"
				[ ${SMPL_PON_SEX_SELECT} == "all_options" ] && PON_SEX1="mixed" && PON_SEX2="M" && PON_SEX3="F" && PON_SEX4="matched_main" && PON_SEX5="matched_each"
				# [ ${SMPL_PON_SEX_SELECT} == "all_options" ] && PON_SEX1="M" && PON_SEX2="mixed" && PON_SEX3="NA" && PON_SEX4="matched_main" && PON_SEX5="matched_each"

				[ ${CTRL} == "no" ] && PON_SEX1="none" && PON_SEX2="NA" && PON_SEX3="NA" && PON_SEX4="NA " && PON_SEX5="NA"

				#LOOP - PoN SEX
				for PON_SEX in ${PON_SEX1} ${PON_SEX2} ${PON_SEX3} ${PON_SEX4} ${PON_SEX5}; do

					#CHECK - PoN SEX
					if [ ${PON_SEX} != "NA" ]; then

						[ ${PON_SEX} == "M" ] && PON_SEX_PATTERN="M" && PON_SEX_ID="M"
						[ ${PON_SEX} == "F" ] && PON_SEX_PATTERN="F" && PON_SEX_ID="F"
						[ ${PON_SEX} == "mixed" ] && PON_SEX_PATTERN="*" && PON_SEX_ID="mixed"
						[ ${PON_SEX} == "matched_main" ] && PON_SEX_PATTERN=${SAMPLE_1_SEX} && PON_SEX_ID=${SAMPLE_1_SEX}
						[ ${PON_SEX} == "matched_each" ] && PON_SEX_PATTERN=${CASE_SEX} && PON_SEX_ID=${CASE_SEX}
						[ ${PON_SEX} == "none" ] && PON_SEX_PATTERN="none" && PON_SEX_ID="none"

						#CHECK - CTRL YES 
						if [ ${CTRL} == "yes" ]; then 
							hdflist=`ls ${CTRL_COUNT_HDF_DIR}*${PON_SEX_PATTERN}"_counts.hdf5" | grep -v ${MAIN_ID} | sed 's/^/-I /'`
							NRealx=`ls ${CTRL_COUNT_HDF_DIR}*${PON_SEX_PATTERN}"_counts.hdf5" | grep -v ${MAIN_ID} | wc -w`
							NReal=$(echo ${NRealx} | tr -d ' ')
						#CHECK - CTRL YES 
						fi

						#CHECK - CTRL YES 
						if [ ${CTRL} == "no" ]; then 
							NReal=0
						#CHECK - CTRL YES 
						fi

						#create PoN
						#============================================================================
						#CHECK - CTRL YES 
						if [ ${CTRL} == "yes" ]; then 
							#CHECK - PON
							if ! ls ${SAMPLE_PoN_DIR}${MAIN_ID}"_PoN-"${PON_SEX_ID}".hdf5" 1> /dev/null 2>&1; then
								echo -e ${COL1}${BEG}${COL2}
								echo -e ${COL1}"◼︎ $(date)"${COL2}
								echo -e ${COL1}"◼︎ "${MAIN_ID}": preparing PoN-"${PON_SEX_ID}${COL2}
								echo -e ${COL1}${END}${COL2}
								echo ${hdflist} > ${SAMPLE_PoN_DIR}${MAIN_ID}"_PoN-"${PON_SEX_ID}".txt"
								${GATK} CreateReadCountPanelOfNormals \
								   ${hdflist} \
								   --annotated-intervals ${INTERVALS_ANNOTATED} \
								   -O ${SAMPLE_PoN_DIR}${MAIN_ID}"_PoN-"${PON_SEX_ID}".hdf5"
							#CHECK - PON
							fi
						#CHECK - CTRL YES 
						fi
						#============================================================================

						#standardization and denoising using %CG correction and PCA by PoN
						#============================================================================
						#CHECK - CTRL YES 
						if [ ${CTRL} == "yes" ]; then 
							# CHECK - PROCESSED DENOIS 
							if ! ls ${RESULTS_DATA_PATTERN}"_denoisedCR_PoN-"${PON_SEX_ID}".tsv" 1> /dev/null 2>&1; then 
								echo -e ${COL1}${BEG}${COL2}
								echo -e ${COL1}"◼︎ $(date)"${COL2}
								echo -e ${COL1}"◼︎ "${MAIN_ID}" - "${CASE_TYPE}" - "${CASE_ID}": denoising read count using PoN-"${PON_SEX_ID}${COL2} 
								echo -e ${COL1}${END}${COL2}	
								${GATK} DenoiseReadCounts \
								    -I ${RESULTS_DATA_PATTERN}"_counts.hdf5" \
								    --count-panel-of-normals ${SAMPLE_PoN_DIR}${MAIN_ID}"_PoN-"${PON_SEX_ID}".hdf5" \
								    --standardized-copy-ratios ${RESULTS_DATA_PATTERN}"_standardizedCR_PoN-"${PON_SEX_ID}".tsv" \
								    --denoised-copy-ratios ${RESULTS_DATA_PATTERN}"_denoisedCR_PoN-"${PON_SEX_ID}".tsv"
							#CHECK - PROCESSED DENOIS 
							fi
						#CHECK - CTRL YES 
						fi 
						#============================================================================


						#standardization and denoising using %CG correction
						#============================================================================
						#CHECK - CTRL YES 
						if [ ${CTRL} == "no" ]; then 
							# CHECK - PROCESSED DENOIS 
							if ! ls ${RESULTS_DATA_PATTERN}"_denoisedCR_PoN-"${PON_SEX_ID}".tsv" 1> /dev/null 2>&1; then 
								echo -e ${COL1}${BEG}${COL2}
								echo -e ${COL1}"◼︎ $(date)"${COL2}
								echo -e ${COL1}"◼︎ "${MAIN_ID}" - "${CASE_TYPE}" - "${CASE_ID}": denoising read count using PoN-"${PON_SEX_ID}${COL2} 
								echo -e ${COL1}${END}${COL2}	
								${GATK} DenoiseReadCounts \
								    -I ${RESULTS_DATA_PATTERN}"_counts.hdf5" \
								    --annotated-intervals ${INTERVALS_ANNOTATED} \
								    --standardized-copy-ratios ${RESULTS_DATA_PATTERN}"_standardizedCR_PoN-"${PON_SEX_ID}".tsv" \
								    --denoised-copy-ratios ${RESULTS_DATA_PATTERN}"_denoisedCR_PoN-"${PON_SEX_ID}".tsv"
							#CHECK - PROCESSED DENOIS 
							fi
						#CHECK - CTRL YES 
						fi 
						#============================================================================

						#segmentation of CN and LOH data
						#============================================================================
						#CHECK - SEGMENTATION 
						if ! ls ${RESULTS_DATA_PATTERN}"_SEGMENTS_PoN-"${PON_SEX_ID}"_"${SEGMENTATIONCONDITIONS}"_segmentation.tsv" 1> /dev/null 2>&1; then
							if ! ls ${RESULTS_DATA_PATTERN}"_SEGMENTS_PoN-"${PON_SEX_ID}"_"${SEGMENTATIONCONDITIONS}".tsv" 1> /dev/null 2>&1; then 
								echo -e ${COL1}${BEG}${COL2}
								echo -e ${COL1}"◼︎ $(date)"${COL2}
								echo -e ${COL1}"◼︎ "${MAIN_ID}" - "${CASE_TYPE}" - "${CASE_ID}": segmentation and quality control, using PoN-"${PON_SEX_ID}${COL2} 
								echo -e ${COL1}${END}${COL2}
								mkdir -p ${QC_DIR}
								
								#generate segmentation contition if smart way is required
								if [[ ${SEGMENTATION_ID} == "smart-"* ]]; then	
									SMART_SEGMENTATION=`Rscript --vanilla ${R}"segmentation_conditions_smart.R" \
									${SEGMENTATION_ID} \
									${RESULTS_DATA_PATTERN}"_denoisedCR_PoN-"${PON_SEX_ID}".tsv" \
									${INTERVALS_ANNOTATED} \
									${SEQ_TYPE}`

									#get rid off the R format
									SMART_SEGMENTATION=$(echo ${SMART_SEGMENTATION} | awk '{print $2}')
									SMART_SEGMENTATION=$(echo ${SMART_SEGMENTATION} | tr -d '\"')
									#extract values that depends on denoised data
									SEGM_DIFFERENCE=$(echo ${SMART_SEGMENTATION} | awk -F"_" '{print $1}') 
				
									SEGMENTATION_FULL=$SEGM_DIFFERENCE$x$SEGM_MINSIZE$x$SEGM_MINSIZE_SUB$x$SEGM_MINPROBE$x$SEGM_MINPROBE_SUB$x$SEGM_MINKEEP$x$SEGM_MINLOSS$x$SEGM_MINLOSS_SUB$x$SEGM_MINLOSS_BIAL$x$SEGM_MINGAIN$x$SEGM_MINGAIN_SUB$x$SEGM_GAP$x$SEGM_SMOOTHPERC$x$SEGM_AFDIF$x$SEGM_AFSIZE$x$SEGM_AFPROBE
								fi

								#segmentation
								Rscript --vanilla ${R}"segmentation_and_qc.R" \
									${UCSC_CHROMOSOME_BANDS} \
									${CAPTURE_ID} \
									${PROJECT_ID} \
									${RESULTS_DATA_PATTERN}"_counts.tsv" \
									${RESULTS_DATA_PATTERN}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".rds" \
									${RESULTS_DATA_PATTERN}"_standardizedCR_PoN-"${PON_SEX_ID}".tsv" \
									${RESULTS_DATA_PATTERN}"_denoisedCR_PoN-"${PON_SEX_ID}".tsv" \
									${MAIN_ID} \
									${CASE_ID} \
									${CASE_TYPE} \
									${CASE_SEX} \
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
									${QCTABLE} \
									${RESULTS_DATA_PATTERN}"_SEGMENTS_PoN-"${PON_SEX_ID}"_"${SEGMENTATIONCONDITIONS}"_segmentation.tsv"
							fi
						#CHECK - SEGMENTATION
						fi
						#============================================================================

						#incorporate gene annotation, ctrl and databases annotation
						#============================================================================
						#CHECK - ANNOTATION DTB
						if ! ls ${RESULTS_DATA_PATTERN}"_SEGMENTS_PoN-"${PON_SEX_ID}"_"${SEGMENTATIONCONDITIONS}".tsv" 1> /dev/null 2>&1; then 
							echo -e ${COL1}${BEG}${COL2}
							echo -e ${COL1}"◼︎ $(date)"${COL2}
							echo -e ${COL1}"◼︎ "${MAIN_ID}" - "${CASE_TYPE}" - "${CASE_ID}": incorporation of gene annotation, DGV and CTRLs to segmentation, using PoN-"${PON_SEX_ID}${COL2} 
							echo -e ${COL1}${END}${COL2}
							#Incorporating databases as DGV, gnomadSV and controls
							Rscript --vanilla ${R}"annotation.R" \
								${DGV} \
								${GENE_ANNOTATION_PROCESSED} \
								${MAPPABILITY} \
								${RESULTS_DATA_PATTERN}"_SEGMENTS_PoN-"${PON_SEX_ID}"_"${SEGMENTATIONCONDITIONS}"_segmentation.tsv" \
								${CTRL_CN_FILE} \
								${N} \
								${CTRL_PON_SEX_SELECT} \
								${SEGM_MINGAIN} \
								${SEGM_MINGAIN_SUB} \
								${SEGM_MINLOSS} \
								${SEGM_MINLOSS_SUB} \
								${RESULTS_DATA_PATTERN}"_SEGMENTS_PoN-"${PON_SEX_ID}"_"${SEGMENTATIONCONDITIONS}".tsv"
						#CHECK - ANNOTATION DTB
						fi
						#============================================================================


					#CHECK - PoN SEX
				    fi

				#LOOP - PoN SEX
				done
				#============================================================================

			#CHECK - SAMPLES EXISTS
			fi

		#====================================================================================
		#LOOP - SAMPLES PROCESSING
		done



		#create family report and plots
		#====================================================================================
		[ ${SMPL_PON_SEX_SELECT} == "M" ] && PON_SEX1="M" && PON_SEX2="NA" && PON_SEX3="NA" && PON_SEX4="NA " && PON_SEX5="NA"
		[ ${SMPL_PON_SEX_SELECT} == "F" ] && PON_SEX1="F" && PON_SEX2="NA" && PON_SEX3="NA" && PON_SEX4="NA" && PON_SEX5="NA"
		[ ${SMPL_PON_SEX_SELECT} == "mixed" ] && PON_SEX1="mixed" && PON_SEX2="NA" && PON_SEX3="NA" && PON_SEX4="NA" && PON_SEX5="NA"
		[ ${SMPL_PON_SEX_SELECT} == "both_separated" ] && PON_SEX1="M" && PON_SEX2="F" && PON_SEX3="NA" && PON_SEX4="NA" && PON_SEX5="NA"
		[ ${SMPL_PON_SEX_SELECT} == "matched_main" ] && PON_SEX1="matched_main" && PON_SEX2="NA" && PON_SEX3="NA" && PON_SEX4="NA" && PON_SEX5="NA"
		[ ${SMPL_PON_SEX_SELECT} == "matched_each" ] && PON_SEX1="matched_each" && PON_SEX2="NA" && PON_SEX3="NA" && PON_SEX4="NA" && PON_SEX5="NA"
		[ ${SMPL_PON_SEX_SELECT} == "all_options" ] && PON_SEX1="mixed" && PON_SEX2="M" && PON_SEX3="F" && PON_SEX4="matched_main" && PON_SEX5="matched_each"
		# [ ${SMPL_PON_SEX_SELECT} == "all_options" ] && PON_SEX1="M" && PON_SEX2="mixed" && PON_SEX3="NA" && PON_SEX4="matched_main" && PON_SEX5="matched_each"
		
		[ ${CTRL} == "no" ] && PON_SEX1="none" && PON_SEX2="NA" && PON_SEX3="NA" && PON_SEX4="NA " && PON_SEX5="NA"

		#LOOP - PoN SEX
		for PON_SEX in ${PON_SEX1} ${PON_SEX2} ${PON_SEX3} ${PON_SEX4} ${PON_SEX5}; do

			#CHECK - PoN SEX
			if [ ${PON_SEX} != "NA" ]; then

				[ ${PON_SEX} == "M" ] && PON_SEX_ID_SAMPLE1="M" && PON_SEX_ID_SAMPLE2="M" && PON_SEX_ID_SAMPLE3="M"
				[ ${PON_SEX} == "F" ] && PON_SEX_ID_SAMPLE1="F" && PON_SEX_ID_SAMPLE2="F" && PON_SEX_ID_SAMPLE3="F"
				[ ${PON_SEX} == "mixed" ] && PON_SEX_ID_SAMPLE1="mixed" && PON_SEX_ID_SAMPLE2="mixed" && PON_SEX_ID_SAMPLE3="mixed"
				[ ${PON_SEX} == "matched_main" ] && PON_SEX_ID_SAMPLE1="${SAMPLE_1_SEX}" && PON_SEX_ID_SAMPLE2="${SAMPLE_1_SEX}" && PON_SEX_ID_SAMPLE3="${SAMPLE_1_SEX}"
				[ ${PON_SEX} == "matched_each" ] && PON_SEX_ID_SAMPLE1="${SAMPLE_1_SEX}" && PON_SEX_ID_SAMPLE2="${SAMPLE_2_SEX}" && PON_SEX_ID_SAMPLE3="${SAMPLE_3_SEX}"
				[ ${PON_SEX} == "none" ] && PON_SEX_ID_SAMPLE1="none" && PON_SEX_ID_SAMPLE2="none" && PON_SEX_ID_SAMPLE3="none"

			 	#create report for each group of samples
				#============================================================================
				#CHECK - REPORT TYPE
				if [ ${PROJECT_TYPE} != "germline" ] && [ ${PROJECT_TYPE} != "tumor" ] && [ ${PROJECT_TYPE} != "compare_independent" ]; then
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ cannot recognized type of the project to generate appropriate report"${COL2}
					echo -e ${COL1}${END}${COL2}
					exit
				fi
				#CHECK - REPORT TYPE

				#CHECK - REPORT
				if [ ${PROJECT_TYPE} == "germline" ] || [ ${PROJECT_TYPE} == "tumor" ]; then
					if ! ls ${RESULTS_DATA_DIR}/${SAMPLE_1_TYPE}_${SAMPLE_1_ID}_${SAMPLE_1_SEX}"_SEGMENTS_PoN-"${PON_SEX}"_"${SEGMENTATIONCONDITIONS}"_REPORT.tsv" 1> /dev/null 2>&1; then 
						echo -e ${COL1}${BEG}${COL2}
						echo -e ${COL1}"◼︎ $(date)"${COL2}
						echo -e ${COL1}"◼︎ "${MAIN_ID}" - "${SAMPLE_1_TYPE}" - "${SAMPLE_1_ID}": creating report, using PoN-"${PON_SEX}${COL2} 
						echo -e ${COL1}${END}${COL2}
						Rscript --vanilla ${R}"create_report.R" \
							${PROJECT_TYPE} \
							${SAMPLE_1_TYPE}"_"${SAMPLE_1_ID} \
							${SAMPLE_1_SEX} \
							${PON_SEX_ID_SAMPLE1} \
							${RESULTS_DATA_DIR}${SAMPLE_1_TYPE}"_"${SAMPLE_1_ID}"_"${SAMPLE_1_SEX}"_SEGMENTS_PoN-"${PON_SEX_ID_SAMPLE1}"_"${SEGMENTATIONCONDITIONS}".tsv" \
							${RESULTS_DATA_DIR}${SAMPLE_1_TYPE}"_"${SAMPLE_1_ID}"_"${SAMPLE_1_SEX}"_denoisedCR_PoN-"${PON_SEX_ID_SAMPLE1}".tsv" \
							${SAMPLE_2_TYPE}"_"${SAMPLE_2_ID} \
							${SAMPLE_2_SEX} \
							${PON_SEX_ID_SAMPLE2} \
							${RESULTS_DATA_DIR}${SAMPLE_2_TYPE}"_"${SAMPLE_2_ID}"_"${SAMPLE_2_SEX}"_denoisedCR_PoN-"${PON_SEX_ID_SAMPLE2}".tsv" \
							${SAMPLE_3_TYPE}"_"${SAMPLE_3_ID} \
							${SAMPLE_3_SEX} \
							${PON_SEX_ID_SAMPLE3} \
							${RESULTS_DATA_DIR}${SAMPLE_3_TYPE}"_"${SAMPLE_3_ID}"_"${SAMPLE_3_SEX}"_denoisedCR_PoN-"${PON_SEX_ID_SAMPLE3}".tsv" \
							${SEGM_MINGAIN} \
							${SEGM_MINGAIN_SUB} \
							${SEGM_MINLOSS} \
							${SEGM_MINLOSS_SUB} \
							${RESULTS_DATA_DIR}${SAMPLE_1_TYPE}_${SAMPLE_1_ID}_${SAMPLE_1_SEX}"_SEGMENTS_PoN-"${PON_SEX}"_"${SEGMENTATIONCONDITIONS}"_REPORT.tsv"
					#CHECK - REPORT
					fi
				fi
				#============================================================================

			 	#create plots for each group of samples
				#============================================================================
				#CHECK - PLOTS
				if ! ls ${RESULTS_PLOT_DIR}"/PoN-"${PON_SEX}"/"${MAIN_ID}"_genome_profile_"${SAMPLE_1_TYPE}"_"${SAMPLE_1_ID}"_"${SAMPLE_1_SEX}".png" 1> /dev/null 2>&1; then
					mkdir -p ${RESULTS_PLOT_DIR}"/PoN-"${PON_SEX}"/detail/"
					mkdir -p ${RESULTS_PLOT_DIR}"/PoN-"${PON_SEX}"/LOH/"
					mkdir -p ${RESULTS_IGV_BED_DIR}"PoN-"${PON_SEX}"/"
					if ls ${R}"plot_CN_"${PROJECT_ID}".R" 1> /dev/null 2>&1; then
						mkdir -p ${RESULTS_PLOT_DIR}"/PoN-"${PON_SEX}"/project_specific_regions/"
					fi
					#Create plots for genome and all chromosomes
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ $MAIN_ID - $SAMPLE_1_TYPE $SAMPLE_1_ID: Creating plots - using PoN-"${PON_SEX}${COL2} 
					echo -e ${COL1}${END}${COL2}
					Rscript --vanilla ${R}"plot_CN.R" \
						${CAPTURE_FILE} \
						${PLOT_ASSEMBLY} \
						${CONTIGS_SIZES} \
						${UCSC_CHROMOSOME_CENTROMERE} \
						${DGV} \
						${MAPPABILITY} \
						${GENE_ANNOTATION_PROCESSED} \
						${PROJECT_ID} \
						${RESULTS_DATA_DIR} \
						${GNOMAD_COUNT_PATTERN} \
						${MAIN_ID} \
						${SAMPLE_1_ID} \
						${SAMPLE_1_TYPE} \
						${SAMPLE_1_SEX} \
						${SAMPLE_2_ID} \
						${SAMPLE_2_TYPE} \
						${SAMPLE_2_SEX} \
						${SAMPLE_3_ID} \
						${SAMPLE_3_TYPE} \
						${SAMPLE_3_SEX} \
						${PON_SEX} \
						${RESULTS_PLOT_DIR}"PoN-"${PON_SEX}"/" \
						${RESULTS_IGV_BED_DIR}"PoN-"${PON_SEX}"/" \
						${SEGMENTATIONCONDITIONS} \
						${SEGM_MINLOSS} \
						${SEGM_MINGAIN} \
				        ${CTRL_CN_FILE} \
				        ${R}"plot_CN_"${PROJECT_ID}".R"
					if ls ${RESULTS_PLOT_DIR}"/PoN-"${PON_SEX}"/detail/"*HOM* 1> /dev/null 2>&1; then
				    	mv ${RESULTS_PLOT_DIR}"/PoN-"${PON_SEX}"/detail/"*HOM* ${RESULTS_PLOT_DIR}"/PoN-"${PON_SEX}"/LOH/"
					fi
					if ls ${RESULTS_PLOT_DIR}"/PoN-"${PON_SEX}"/detail/"*OTHER* 1> /dev/null 2>&1; then
				    	mv ${RESULTS_PLOT_DIR}"/PoN-"${PON_SEX}"/detail/"*OTHER* ${RESULTS_PLOT_DIR}"/PoN-"${PON_SEX}"/LOH/"
					fi
					#transform bedGraph to BigWig
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ $MAIN_ID: transforming bedGraph to BigWig - using PoN "${PON_SEX}${COL2} 
					echo -e ${COL1}${END}${COL2}
					
					for file in ${RESULTS_IGV_BED_DIR}"PoN-"${PON_SEX}/*"_SNP"; do
						${BEDGRAPHTOBIGWIG} ${file} ${CONTIGS_SIZES} ${file}".bw"
						rm ${file}
					done
					for file in ${RESULTS_IGV_BED_DIR}"PoN-"${PON_SEX}/*"_CN"; do
						${BEDGRAPHTOBIGWIG} ${file} ${CONTIGS_SIZES} ${file}".bw"
						rm ${file}
					done
				
				#CHECK - PLOTS
				fi
				#============================================================================
			

			 	#create plots for each group of samples - detail for regions of interest
				#============================================================================
				REGIONS_TABLE=`cat ${MASTERREGOFINT} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" -v MAIN_ID="${MAIN_ID}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION && $5==MAIN_ID)print}'`

				if [ -n "${REGIONS_TABLE}" ]; then

					#remove term table if exists
					if ls ${MASTERREGOFINT}"_temp.txt" 1> /dev/null 2>&1; then
						rm ${MASTERREGOFINT}"_temp.txt"
					fi
					#create header
					echo -e "CONTIG\tSTART\tEND\tZOOM" > ${MASTERREGOFINT}"_temp.txt"
					#filter for regions of interest
					cat ${MASTERREGOFINT} | awk -F"\t" -v PROJECT_ID="${PROJECT_ID}" -v CAPTURE_ID="${CAPTURE_ID}" -v GENOME_VERSION="${GENOME_VERSION}" -v MAIN_ID="${MAIN_ID}" '{if ($1=="yes" && $2==PROJECT_ID && $3==CAPTURE_ID && $4==GENOME_VERSION && $5==MAIN_ID)print}' | while read -r line || [[ -n "$line" ]]; do
						CONTIG=`extract_from_master 'CONTIG' ${MASTERREGOFINT}`
						START=`extract_from_master 'START' ${MASTERREGOFINT}`
						END=`extract_from_master 'END' ${MASTERREGOFINT}`
						ZOOM=`extract_from_master 'ZOOM' ${MASTERREGOFINT}`

						echo -e ${CONTIG}'\t'${START}'\t'${END}'\t'${ZOOM} >> ${MASTERREGOFINT}"_temp.txt"
						# echo $line >> ${MASTERREGOFINT}"_temp.txt"
					done

					mkdir -p ${RESULTS_PLOT_DIR}"/PoN-"${PON_SEX}"/detail_regions_of_interest/"
					
					#Create plots for regions of interest
					echo -e ${COL1}${BEG}${COL2}
					echo -e ${COL1}"◼︎ $(date)"${COL2}
					echo -e ${COL1}"◼︎ $MAIN_ID - $SAMPLE_1_TYPE $SAMPLE_1_ID: Creating additional plots for REGIONS of INTEREST - using PoN-"${PON_SEX}${COL2}
					echo -e ${COL1}${END}${COL2}

					Rscript --vanilla ${R}"plot_CN_reg_of_interest.R" \
						${CAPTURE_FILE} \
						${PLOT_ASSEMBLY} \
						${CONTIGS_SIZES} \
						${UCSC_CHROMOSOME_CENTROMERE} \
						${DGV} \
						${MAPPABILITY} \
						${GENE_ANNOTATION_PROCESSED} \
						${PROJECT_ID} \
						${RESULTS_DATA_DIR} \
						${GNOMAD_COUNT_PATTERN} \
						${MAIN_ID} \
						${SAMPLE_1_ID} \
						${SAMPLE_1_TYPE} \
						${SAMPLE_1_SEX} \
						${SAMPLE_2_ID} \
						${SAMPLE_2_TYPE} \
						${SAMPLE_2_SEX} \
						${SAMPLE_3_ID} \
						${SAMPLE_3_TYPE} \
						${SAMPLE_3_SEX} \
						${PON_SEX} \
						${RESULTS_PLOT_DIR}"PoN-"${PON_SEX}"/" \
						${RESULTS_IGV_BED_DIR}"PoN-"${PON_SEX}"/" \
						${SEGMENTATIONCONDITIONS} \
						${SEGM_MINLOSS} \
						${SEGM_MINGAIN} \
				        ${CTRL_CN_FILE} \
				        ${R}"plot_CN_"${PROJECT_ID}".R" \
				        ${MASTERREGOFINT}"_temp.txt"

				    rm ${MASTERREGOFINT}"_temp.txt"

				fi
				#============================================================================

			#CHECK - PoN SEX
		    fi
		
		#====================================================================================
		# LOOP PON FOR REPORT AND PLOTS
		done


	#LOOP - LINE oF MASTER
	done
	#============================================================================

	# remove temp files
	SEGMENTATION_TEMP=${SEGMENTATION}"_temp.txt"
	if ls ${SEGMENTATION_TEMP} 1> /dev/null 2>&1; then
		rm ${SEGMENTATION_TEMP}
	fi
	
#LOOP - LINE oF PROJECTS
done
#============================================================================

echo -e ${COL1}${BEG}${COL2}
echo -e ${COL1}"◼︎ $(date)"${COL2}
echo -e ${COL1}"◼︎ samples analysis finished"${COL2}
echo -e ${COL1}${END}${COL2}

####============================================================================


