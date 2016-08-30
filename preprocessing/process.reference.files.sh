#!/bin/bash

# manipulate gencode files in various ways to make them easier to use downstream
# takes in GTEx and gencode annotation files from the command line

set -o nounset -o errexit -o pipefail

if [ $# -ne 2 ]; then
	echo "Usage: process.reference.files.sh <gtex annotation, gzipped> <gencode annotation>"
	exit
fi

gtex=$1
gencode=$2

gencodeprefix=${RAREVARDIR}/reference/gencode.v19.annotation
gtexprefix=${RAREVARDIR}/reference/gencode.v19.genes.v6p.patched_contigs

# TSS
bash gtf2TSS.sh <(zcat $gtex) > ${gtexprefix}_TSS.bed

# gene bed file
bash gtf2genebed.sh -t <(zcat $gtex) | awk '{print "chr"$0}' > ${gtexprefix}.bed

# genetypes
zcat $gtex | awk '{if($3=="transcript" && $1 ~ /^[0-9]+$/){print substr($10,2,length($10)-3)"\t"substr($14,2,length($14)-3)}}' \
    > ${gtexprefix}_genetypes_autosomal.txt

# for gencode annotation get lincRNA and prootein-coding exons, but also pad the internal exons by 5 base pairs
python pad.gtf.exons.py $gencode | sort -k1,1 -k2,2n | uniq > ${gencodeprefix}_coding.lincRNA_padded.bed
