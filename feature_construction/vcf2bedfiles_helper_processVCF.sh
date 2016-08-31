#!/bin/bash

# uses vcftools to extract relevant information from vcf files
# helper script for vcf2bedfiles

# take as input:
# * vcf file
# * file with individuals to include
# * output directory

if [ $# -ne 3 ]; then
    echo "usage: vcf2bedfiles_helper_provessVCF.sh <vcf> <individuals file> <output prefix>"
    exit
fi

vcf=$1
indincl=$2 # this should be only european american individuals
outprefix=$3

# individuals to include from the allele frequency calculation
# generate file name from the indincl file name
keeplist=${indincl%_ids.txt}_VCFids.txt # as above, but IDs as they appear in the snp/indel VCF

# get their ids as shown in the vcf file
zcat $vcf | awk '{if(substr($1,1,2)!="##"){print; exit}}' | \
awk -v indincl=$indincl 'BEGIN{
    while((getline<indincl)>0){
        IDs[$1]
    }
}{
    for(i=10;i<=NF;i++){
        split($i,IDparts,"-"); 
        indID=IDparts[1]"-"IDparts[2]; 
        if(indID in IDs){print $i}
    }
}' > $keeplist

# SNPs
vcftools --gzvcf $vcf --out ${outprefix}_SNPs --remove-filtered-all --keep $keeplist --remove-indels --max-missing-count 10 --freq &
vcftools --gzvcf $vcf --out ${outprefix}_SNPs --remove-filtered-all --keep $keeplist --remove-indels --max-missing-count 10 --extract-FORMAT-info GT &

# indels
vcftools --gzvcf $vcf --out ${outprefix}_indels --remove-filtered-all --keep $keeplist --keep-only-indels --max-missing-count 10 --freq &
vcftools --gzvcf $vcf --out ${outprefix}_indels --remove-filtered-all --keep $keeplist --keep-only-indels --max-missing-count 10 --extract-FORMAT-info GT

wait
