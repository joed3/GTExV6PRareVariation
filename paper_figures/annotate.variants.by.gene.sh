#!/bin/bash

set -o nounset -o errexit -o pipefail

# Set the path to the variant files
VARS_DIR=${RAREVARDIR}/features/variantBeds/individuals/

# Set the path to the gene annotation
GENE=${RAREVARDIR}/reference/gencode.v19.genes.v6p.patched_contigs.bed

# Set the path to the ER promoter annotation
PROM_DIR=${RAREVARDIR}/features/annotations/epigenomicsRoadmap/prom

# Set output directory
OUT=${RAREVARDIR}/features/variantBeds/gtex.vars.maf0-25.er.prom.bed

#-- Process variant data
# Read in variant data for each individual
# Collapse across individuals
# Overlap with gene data
# Overlap with promoter data from ER 

if [ -e ${OUT} ]; then
	rm ${OUT}
fi 

for vartype in "SNPs" "indels" "HallLabSV"; do 
	echo ${vartype}
	if [ ${vartype} == "HallLabSV" ]; then
		distance=200000
	else
		distance=10000
	fi
	echo ${distance}
	zcat ${VARS_DIR}/*_${vartype}.bed.gz | cut -f1,2,3,4 | sort -k1,1 -k2,2n | uniq | \
		bedtools window -r ${distance} -l ${distance} -a ${GENE} -b stdin | \
		awk -v vartype=${vartype} 'BEGIN{OFS="\t"}{print $5, $6, $7, $8, $4, vartype}' | \
		sort -k1,1 -k2,2n | \
		bedtools intersect -sorted -wa -c -a stdin -b ${PROM_DIR}/*.bed.gz >> ${OUT}
done	

# Combine the variant files
bgzip ${OUT}

