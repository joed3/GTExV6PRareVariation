#!/bin/bash

set -o nounset -o errexit -o pipefail

# set number of processes
nproc=10

## paths are hardcoded and woul dneed to be updated for a different version/purpose!
indir=${RAREVARDIR}/features/byGene
inds=${RAREVARDIR}/preprocessing/gtex_2015-01-12_wgs_ids_outlier_filtered.txt
outlierdir=${RAREVARDIR}/data/singlez

# set tissue to use
tissues=`cut -f1 ${RAREVARDIR}/preprocessing/gtex_2015-01-12_tissue_by_ind.txt | tail -n+2 | uniq | tr '\n' ' '`

# function takes as input:
# * the window over which the variants were aggregated
# * the tissue
runfeatures() {
    window=$1
    tissue=$2
    echo "counts $window $tissue"
    python pick_outliers_controls_imbalanced.py \
	--FEATURE_DIR ${indir}/${window}/counts/MAF0-1/SNPs \
	--OUTLIER_PICKED ${outlierdir}/outliers_singlez_nothreshold_${tissue}_picked.txt \
	--COUNTS ${outlierdir}/outliers_singlez_nothreshold_${tissue}_counts.txt \
	--INDS $inds \
	--OUT ${outlierdir}/${window}_features_SNPs_MAF0-1_nothreshold_${tissue}.txt \
	--type Z \
	--threshold 0 \
	--count_threshold 1 \
	--END "_counts_bygene.txt"
}

# make all the necessary variables avaialble to parallel
export -f runfeatures
export indir
export outlierdir
export inds

parallel --jobs $nproc runfeatures ::: 10kb ::: $tissues

wait
echo "done compiling features for singlez"

