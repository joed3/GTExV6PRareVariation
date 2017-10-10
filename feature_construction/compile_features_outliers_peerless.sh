#!/bin/bash

# author: Emily Tsang
# modified by Joe Davis

set -o nounset -o errexit -o pipefail

## takes as input the suffix indicating which PEER corrected data to use
## assumes that want features for all subdirectories of the input directory with the exception of those that start with archive
## also assumes those directories will have the correct structure (i.e. will have been created by run_build_feature_count_summaries_all_genes.sh)
## will make a directory based on the method name does not exist.
                                                     
usage="usage: run_compile_features_outliers_peerless.sh <PEERSuffix>"
if [ $# -ne 1 ]; then
    echo "Error: Must supply the suffix indicating which PEER-corrected data to use."
    echo $usage
    exit
fi

PEERSuffix=$1

# Define global variables
METHOD=medz
dir=${RAREVARDIR}/features
inds_feat=${RAREVARDIR}/preprocessing/gtex_2015-01-12_wgs_ids_outlier_filtered${PEERSuffix}
inds_count=${RAREVARDIR}/preprocessing/gtex_2015-01-12_wgs_ids_HallLabSV_outlier_filtered${PEERSuffix}
indir=${dir}/byGene
outdir=${RAREVARDIR}/data
outlierdir=${outdir}/${METHOD}
outliers=${outdir}/outliers_${METHOD}_picked${PEERSuffix}
counts=${outdir}/outliers_${METHOD}_counts${PEERSuffix}
stattype=Z
threshold=2

# function takes as input:
# * the variant type (SNPs, indels or HallLabSV)
# * the MAF bin
# * the window over which the variants were aggregated
# * "T" if we are dealing with count features
runfeatures() {
    vartype=$1
    maf=$2
    window=$3
    is_counts=$4
    if [ $is_counts == "T" ]; then
	inputdir=${indir}/${window}/counts/MAF${maf}/${vartype}
	inds=$inds_count
	output_file=${outlierdir}/${window}_features_${vartype}_counts_MAF${maf}${PEERSuffix}
	suffix="_counts_bygene.txt"
	echo "counts $window $vartype $maf"
    else
	inputdir=${indir}/${window}/MAF${maf}/${vartype}
	inds=$inds_feat
	output_file=${outlierdir}/${window}_features_${vartype}_MAF${maf}${PEERSuffix}
	suffix="_features_bygene.txt"
	echo "features $window $vartype $maf"
    fi
    python feature_construction/pick_outliers_controls_imbalanced.py \
	--FEATURE_DIR $inputdir \
	--OUTLIER_PICKED $outliers \
	--COUNTS $counts \
	--INDS $inds \
	--OUT $output_file \
	--type $stattype \
	--threshold $threshold \
	--END $suffix
}

# make all the necessary variables avaialble to parallel
export -f runfeatures
export PEERSuffix
export indir
export outlierdir
export inds_count
export inds_feat
export stattype
export threshold
export outliers
export counts

# make the features for all window sizes (each a subdirectory of in the input directory)
subdirs=`ls $indir`
windows=""
for subdir in $subdirs
do
    # skip directory names starting with archive
    if [[ $subdir != "archive"* ]]; then
	windows="$subdir $windows"
    fi
done

# annotation features
parallel --jobs 10 runfeatures ::: SNPs indels ::: 0-1 1-5 5-10 10-25 ::: $windows ::: F
# count features
parallel --jobs 10 runfeatures ::: SNPs indels HallLabSV ::: 0-1 1-5 5-10 10-25 ::: $windows ::: T


echo "done compiling features for $METHOD!"

