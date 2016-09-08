#!/bin/bash

# compile the site-level feature for each genes
# only include features for individuals that have at least one variant within the vicinity of the gene

set -o nounset -o errexit -o pipefail

# takes as input several parameters:
# * base directory in which to write; will create this directory with subdirectories
# * window size
# * bed file with intervals around which to build the window (either with TSSs or entire gene bodies)
# * bed file with coding regions to exclude (optional)

### SET THE NUMBER OF PROCESSES TO RUN IN PARALLEL
nproc=20

if [ $# -ge 4 ] && [ $# -le 5 ]; then
    IN_DIR=$1
    OUT_DIR=$2
    WINDOW_SIZE=$3
    GENE=$4
    if [ $# -eq 5 ]; then
	CODING=$5
    else
	# make dummy file that is empty
	if [ ! -e dummy.bed ]; then
	    echo -e "chr1\t0\t1" > dummy.bed
	fi
	CODING=dummy.bed
    fi
else
    echo "usage: build_feature_summaries_all_genes.sh site_features_bed.gz output_directory window_size gene_bed (exclude_bed)"
    exit
fi

# make output directory (and subdirectories) if not already there
# if the upper-level directory is there, assume the subdirectories are also there
current_dir=`pwd`
if [ ! -d $OUT_DIR ]; then
    mkdir $OUT_DIR
    cd $OUT_DIR
    mkdir MAF0-1
    mkdir MAF1-5
    mkdir MAF5-10
    mkdir MAF10-25
    for dir in MAF*
    do
	cd $dir
	mkdir SNPs
	mkdir indels
	cd ..
    done
fi
cd $current_dir

# get window of variants around each gene and feed that into feature building script
# do this for different MAF thresholds
# MAF 0-1%, 1-5%, 5-10%, 10-25% where the beginning of the range is excluded and the end is included
makegenefeatures() {
    var_features=$1

    endpoints=(0 0.01 0.05 0.1 0.25)
    names=(0 1 5 10 25)

    # strip features bed file name to get gtex id and variant type
    fname=`basename $var_features`
    sample=${fname%%_*}
    vartype=${fname%_features.bed.gz}
    vartype=${vartype#*_}
    
    # create output file name
    out_prefix=${sample}_features_bygene
    
    for ((i=1; i<${#endpoints[@]}; i++))
    do
	mafname=${names[$i-1]}-${names[$i]}
	zcat $var_features | \
	    awk -v start=${endpoints[$i-1]} -v end=${endpoints[$i]} '$4>start && $4<=end && (end>0.01 || $(NF-1)=="NA" || $(NF-1)<=0.01)' | \
	    bedtools intersect -sorted -wa -v -a stdin -b $CODING | \
	    bedtools window -r $WINDOW_SIZE -l $WINDOW_SIZE -a $GENE -b stdin | \
	    python build_feature_set.py --features ${OUT_DIR}/MAF${mafname}/${vartype}/${out_prefix}.txt
    done
}

# make things available to parallel (note: cannot export arrays, so endpoints and names are defined in the function)
export -f makegenefeatures
export CODING
export WINDOW_SIZE
export GENE
export OUT_DIR

parallel --jobs $nproc makegenefeatures ::: ${IN_DIR}/GTEX-*

echo "all done!"
