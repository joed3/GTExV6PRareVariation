#!/bin/bash

## Pipeline to build TFBS BED files for feature creation
## Will be using the TFBSs from Kheradpour and Kellis 2013


## Make directory in the goats annotations directory to hold TFBS data
outdir=${RAREVARDIR}/features/annotations/TFBS_Pouya

if [ ! -d $outdir ]; then
    mkdir $outdir
fi

## Create BED file of raw TFBSs
zcat ${ENCODE_MOTIF_DIR}/matches.txt.gz | grep 'known\|disc' | sed 's/ /\t/g' \
	| awk '{OFS="\t"}{split($1, a, "_"); print $2, $3-1, $4+1, $1, $5, a[1]}' | \
	gzip - > ${outdir}/matches.bed.gz

## Get summary stats like number of motifs per TF, length, number of TFBSs, etc.
## Make summary plots
python pouya.raw.summary.py ${outdir}/matches.bed.gz ${outdir}/matches.raw.summary.txt

## Splits the processed motifs file from the Pouya dataset into BED files for each transcription factor
tail -n +2 ${outdir}/matches.raw.summary.txt | cut -f1 | sort | uniq > ${outdir}/tfs.txt
for line in $(cat ${outdir}/tfs.txt)
do
    echo ${line}
    zcat ${outdir}/matches.bed.gz | grep ${line} \
    | grep -v 'chrX\|chrY' | sort -k1,1 -k2,2n \
    | gzip - > ${outdir}/${line}.bed.gz
done

## Merges overlapping regions for each of the TF BED files
for line in $(cat ${outdir}/tfs.txt)
do
    echo ${line}
    bedtools merge -i ${outdir}/${line}.bed.gz \
    | gzip - > ${outdir}/${line}.merged.bed.gz
done