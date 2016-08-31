#!/bin/bash

set -o nounset -o errexit -o pipefail

## Required file: er.tissues.map.grouped.tsv

## Copies selected ER tissues to our annotation directory
## Also makes sure the files are position sorted
ANNO_DIR=${RAREVARDIR}/features/annotations/epigenomicsRoadmap

if [ ! -d $ANNO_DIR ]; then
    mkdir $ANNO_DIR
fi

for anno in dyadic enh prom
do
	if [ ! -d ${ANNO_DIR}/${anno} ]; then
		mkdir ${ANNO_DIR}/${anno}
	fi
	header=1
	while IFS='' read -r line || [[ -n "$line" ]]
	do
	    if [ $header -eq 0 ]
	    then
		arrLine=( $line )
		ID=${arrLine[0]}
		TISSUE=${arrLine[1]}
		zcat ${ER_DIR}/${anno}/BED_files_per_sample/regions_${anno}_${ID}.bed.gz | sort -k1,1 -k2,2n | \
		gzip - > ${ANNO_DIR}/${anno}/${TISSUE}.bed.gz 
	    else
		header=0
	    fi
	done < ${RAREVARDIR}/preprocessing/er.tissues.map.grouped.tsv
done

## Merges similar tissues into single BED file
## Collapses overlapping features

## Brain
for anno in dyadic enh prom
do
	zcat ${ANNO_DIR}/${anno}/Brain_Angular_Gyrus.bed.gz ${ANNO_DIR}/${anno}/Brain_Anterior_Caudate.bed.gz \
	${ANNO_DIR}/${anno}/Brain_Cingulate_Gyrus.bed.gz ${ANNO_DIR}/${anno}/Brain_Hippocampus_Middle.bed.gz \
	${ANNO_DIR}/${anno}/Brain_Inferior_Temporal_Lobe.bed.gz ${ANNO_DIR}/${anno}/Brain_Mid_Frontal_Lobe.bed.gz \
	| sort -k1,1 -k2,2n | bedtools merge -i - | gzip - > ${ANNO_DIR}/${anno}/Brain.bed.gz
	rm ${ANNO_DIR}/${anno}/Brain_Angular_Gyrus.bed.gz ${ANNO_DIR}/${anno}/Brain_Anterior_Caudate.bed.gz \
	${ANNO_DIR}/${anno}/Brain_Cingulate_Gyrus.bed.gz ${ANNO_DIR}/${anno}/Brain_Hippocampus_Middle.bed.gz \
	${ANNO_DIR}/${anno}/Brain_Inferior_Temporal_Lobe.bed.gz ${ANNO_DIR}/${anno}/Brain_Mid_Frontal_Lobe.bed.gz 
done

## Colon
for anno in dyadic enh prom
do
	zcat ${ANNO_DIR}/${anno}/Colonic_Mucosa.bed.gz ${ANNO_DIR}/${anno}/Colon_Smooth_Muscle.bed.gz \
	| sort -k1,1 -k2,2n | bedtools merge -i - | gzip - > ${ANNO_DIR}/${anno}/Colon_Transverse.bed.gz
	rm ${ANNO_DIR}/${anno}/Colonic_Mucosa.bed.gz ${ANNO_DIR}/${anno}/Colon_Smooth_Muscle.bed.gz
done

## Muscle Skeletal
for anno in dyadic enh prom
do
	zcat ${ANNO_DIR}/${anno}/Skeletal_Muscle_Female.bed.gz ${ANNO_DIR}/${anno}/Skeletal_Muscle_Male.bed.gz \
	| sort -k1,1 -k2,2n | bedtools merge -i - | gzip - > ${ANNO_DIR}/${anno}/Muscle_Skeletal.bed.gz
	rm ${ANNO_DIR}/${anno}/Skeletal_Muscle_Female.bed.gz ${ANNO_DIR}/${anno}/Skeletal_Muscle_Male.bed.gz
done

## Stomach
for anno in dyadic enh prom
do
	zcat ${ANNO_DIR}/${anno}/Stomach_Mucosa.bed.gz ${ANNO_DIR}/${anno}/Stomach_Smooth_Muscle.bed.gz \
	| sort -k1,1 -k2,2n | bedtools merge -i - | gzip - > ${ANNO_DIR}/${anno}/Stomach.bed.gz
	rm ${ANNO_DIR}/${anno}/Stomach_Mucosa.bed.gz ${ANNO_DIR}/${anno}/Stomach_Smooth_Muscle.bed.gz
done
