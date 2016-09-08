#!/bin/bash

# Bash script to pick outliers and controls for MEDZ with threshold = 0
# In other words, picks the most extreme outlier per gene and considers all remaining tested individuals as controls

set -o nounset

# Run python script to pick most extreme individuals per gene as outliers 
# All others as controls
FEATURE_DIR=${RAREVARDIR}/features/byGene
INDS_COUNT=${RAREVARDIR}/preprocessing/gtex_2015-01-12_wgs_ids_outlier_filtered.txt
OUTLIER_DIR=${RAREVARDIR}/data
TYPE=SNPs

python pick_outliers_controls_imbalanced.py --FEATURE_DIR ${FEATURE_DIR}/10kb/counts/MAF0-1/${TYPE} \
	--OUTLIER_PICKED ${OUTLIER_DIR}/outliers_medz_nothreshold_picked.txt \
	--COUNTS ${OUTLIER_DIR}/outliers_medz_counts.txt \
	--INDS ${INDS_COUNT} \
	--OUT ${OUTLIER_DIR}/medz/10kb_features_${TYPE}_counts_MAF0-1_nothreshold.txt --type Z --threshold 0 \
	--END _counts_bygene.txt

