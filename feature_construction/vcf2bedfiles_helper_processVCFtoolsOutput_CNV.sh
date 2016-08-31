#!/bin/bash

set -o nounset
shopt -s expand_aliases

# helper script for vcf2bedfiles.sh
# processes CNV output
# this is done separately from other SVs because the genotypes are encoded differently
# uses genotype frequencies instead of allele frequencies

# NOTE: appends directly to the HallLabSV bedfiles
# therefore, make sure to not run this script multiple times in a row without running the other helper scripts

# encodes all individuals without the major genotype as heterozygous

# required files:
# GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller_123EAonly_HallLabCNV.CN.FORMAT
# GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller_123EAonly_HallLabCNV.INFO

# takes as input the output directory
usage="usage: vcf2bedfiles_helper_processVCFtoolsOutput_CNV.sh <I/O directory>"
if [ $# -ne 1 ]; then
    echo $usage
    exit
fi

# FILE PATHS
dir=$1
cn=${dir}/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller_123EAonly_HallLabCNV.CN.FORMAT
info=${dir}/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller_123EAonly_HallLabCNV.INFO
af=${dir}/AF_HallLabCNV.bed

# SOME FUNCTIONS                                                                                    
# skip first line of file                                                                           
alias skip_header="tail -n+2"

# MAKE ALLELE FREQUENCY FILE
# columns: CHROM POS AF CN (where CN is the most common CN)
# if multiple CN have the same count, pick the smaller one
cat $cn | skip_header | \
awk 'BEGIN{OFS="\t"; print "CHROM","POS","GT_FREQ","MOST_COMMON_CN"}{
  delete arr;
  for(i=3; i<=NF; i++) {arr[$i]=arr[$i]+1};
  maxCN = -1;
  CNcount = 0;
  for(cn in arr) {
    if(arr[cn] > CNcount || (arr[cn]==CNcount && cn < maxCN)) {CNcount=arr[cn]; maxCN=cn} 
  };
  af = 1 - CNcount/(NF-2);
  print $1,$2,af,maxCN;
}' > $af

# PROCESS COPY NUMBER FILE
# use end coordinate information from the INFO file
processcn() {
    i=$1
    sample=`head -n1 $cn | cut -f$i`
    paste <(cut -f1,2,$i $cn) $af $info | skip_header | \
	awk 'BEGIN{OFS="\t"}{
               if ($6>0.25){next};
               if ($1=="X" || $1=="Y"){next};
               if($3!=$7){print "chr"$1,$2-1,$12,$6,1}}' \
    >> ${dir}/individuals/${sample}_HallLabSV.bed
    # sort
    sort -k1,1 -k2,2n ${dir}/individuals/${sample}_HallLabSV.bed > ${dir}/${sample}.tmp
    mv -f ${dir}/${sample}.tmp ${dir}/individuals/${sample}_HallLabSV.bed
}

# make things available to parallel
export -f processcn
export cn
export af
export dir

ncol=`head -n1 $cn | wc -w`
cols=$(eval echo "{3..$ncol}")

for i in $cols; do echo $i; done | parallel --jobs 10 processcn {}

wait # just to be sure it doesn't return before it's done.
