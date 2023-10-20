#!/bin/bash

ABSCN_DATA_DIR=${RESULTS_DIR_GO}"ascat_C-"${N}"_PoN-"${PON_SEX}"_"${SEGMENTATIONCONDITIONS}"/data/"
#ABSCN_PLOT_DIR=${RESULTS_DIR_GO}"ascat_C-"${N}"_PoN-"${PON_SEX}"_"${SEGMENTATIONCONDITIONS}"/plot/"
ABSCN_IGV_DIR=${RESULTS_DIR_GO}"ascat_C-"${N}"_PoN-"${PON_SEX}"_"${SEGMENTATIONCONDITIONS}"/IGV/"

OUTPUT=${ABSCN_DATA_DIR}${SAMPLE_1_TYPE}"_"${SAMPLE_1_ID}"_"${SAMPLE_1_SEX}"_SEGMENTS_C-"${N}"_PoN-"${PON_SEX_ID_SAMPLE1}"_"${SEGMENTATIONCONDITIONS}"_ascat.rds"
#CHECK - ASCAT
if ! ls ${OUTPUT} 1> /dev/null 2>&1; then 

	echo -e ${COL1}${BEG}${COL2}
	echo -e ${COL1}"◼︎ $(date)"${COL2}
	echo -e ${COL1}"◼︎ $MAIN_ID - $SAMPLE_1_TYPE $SAMPLE_1_ID: Running ascat - using PoN-"${PON_SEX}${COL2} 
	echo -e ${COL1}${END}${COL2}
	mkdir -p ${ABSCN_DATA_DIR}
	mkdir -p ${ABSCN_IGV_DIR}
	Rscript --vanilla ${R}"ascat.R" \
		${ASCAT_FRACTION} \
		${SAMPLE_1_TYPE}"_"${SAMPLE_1_ID} \
		${RESULTS_DATA_DIR}${SAMPLE_1_TYPE}"_"${SAMPLE_1_ID}"_"${SAMPLE_1_SEX}"_SEGMENTS_C-"${N}"_PoN-"${PON_SEX_ID_SAMPLE1}"_"${SEGMENTATIONCONDITIONS}"_segmentation.tsv" \
		${RESULTS_DATA_DIR}${SAMPLE_1_TYPE}"_"${SAMPLE_1_ID}"_"${SAMPLE_1_SEX}"_denoisedCR_C-"${N}"_PoN-"${PON_SEX_ID_SAMPLE1}".tsv" \
		${RESULTS_DATA_DIR}${SAMPLE_1_TYPE}"_"${SAMPLE_1_ID}"_"${SAMPLE_1_SEX}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".rds" \
		${RESULTS_DATA_DIR}${SAMPLE_2_TYPE}"_"${SAMPLE_2_ID}"_"${SAMPLE_2_SEX}"_denoisedCR_C-"${N}"_PoN-"${PON_SEX_ID_SAMPLE2}".tsv" \
		${RESULTS_DATA_DIR}${SAMPLE_2_TYPE}"_"${SAMPLE_2_ID}"_"${SAMPLE_2_SEX}"_allelicCounts"${GNOMAD_COUNT_PATTERN}".rds" \
		${ABSCN_DATA_DIR}${SAMPLE_1_TYPE}"_"${SAMPLE_1_ID}"_"${SAMPLE_1_SEX}"_SEGMENTS_C-"${N}"_PoN-"${PON_SEX_ID_SAMPLE1}"_"${SEGMENTATIONCONDITIONS}"_segmentation_absCN.tsv" \
		${ABSCN_DATA_DIR}${SAMPLE_1_TYPE}"_"${SAMPLE_1_ID}"_"${SAMPLE_1_SEX}"_denoisedCR_C-"${N}"_PoN-"${PON_SEX_ID_SAMPLE1}"_absCN.tsv" \
		${OUTPUT} \
		${ABSCN_DATA_DIR} \
		${ABSCN_IGV_DIR}

# PURITY=`echo $VARIABLE | cut -d'_' -f1`
# PLOIDY=`echo $VARIABLE | cut -d'_' -f2`

	#transform bedGraph to BigWig
	echo -e ${COL1}${BEG}${COL2}
	echo -e ${COL1}"◼︎ $(date)"${COL2}
	echo -e ${COL1}"◼︎ $MAIN_ID: ASCAT - transforming bedGraph to BigWig - using PoN "${PON_SEX}${COL2} 
	echo -e ${COL1}${END}${COL2}
	
	# for file in ${ABSCN_IGV_DIR}*"_CN_absCN"; do
	# 	${BEDGRAPHTOBIGWIG} ${file} ${CONTIGS_SIZES} ${file}".bw"
	# 	rm ${file}
	# done
	for file in ${ABSCN_IGV_DIR}*"_segments_absCN"; do
		${BEDGRAPHTOBIGWIG} ${file} ${CONTIGS_SIZES} ${file}".bw"
		rm ${file}
	done

	for file in ${ABSCN_IGV_DIR}*"_allele"; do
		${BEDGRAPHTOBIGWIG} ${file} ${CONTIGS_SIZES} ${file}".bw"
		rm ${file}
	done

fi



