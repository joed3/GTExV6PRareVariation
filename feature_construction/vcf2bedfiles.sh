#!/bin/bash

# Goes from GTEx VCF files to a bed file for each individual with that individual's variant sites and their allele frequency

set -o nounset

# set number of processes to run in parallel when gzipping
nproc=15

## takes as input the two vcfs to operate on
usage="usage: vcf2bedfiles.sh <GTEx SNV/indel vcf> <SV vcf>"
if [ $# -ne 2 ]; then
    echo $usage
    exit
fi

vcf=$1
vcfsv=$2

indincl=../preprocessing/gtex_2015-01-12_AF_ids.txt
beddir=${RAREVARDIR}/features/variantBeds

# get prefix of output files
nEA=`wc -l $indincl | awk '{print $1}'`
fileprefix=`basename $vcf`
fileprefix=${fileprefix%.vcf.gz}
fileprefix=${fileprefix}"_"${nEA}"EAonly"

prefix=${beddir}/${fileprefix}

## actually run things!
#######################
date
# first get vctools to generate useful information
echo "Processing VCF files with vcftools..."
bash vcf2bedfiles_helper_processVCF.sh $vcf $indincl $prefix
echo "Processing VCF files (SNPs/indels) done."
date

# process SNPs
echo
echo "Processing SNPs..."
bash vcf2bedfiles_helper_processVCFtoolsOutput.sh SNPs $prefix
sleep 5 # so they don't both create the outdir at the same time
echo
# process indels
echo "Processing indels..."
bash vcf2bedfiles_helper_processVCFtoolsOutput.sh indels $prefix
date

# also run for SVs
echo
echo "Processing SV VCF file..."
bash vcf2bedfiles_helper_processVCF_SV.sh $vcfsv $beddir
echo "Processing SV VCF file done."
date
echo
echo "Processing SVs..."
bash vcf2bedfiles_helper_processVCFtoolsOutput.sh HallLabSV $prefix
# process CNVs (append to SV output files)
bash vcf2bedfiles_helper_processVCFtoolsOutput_CNV.sh $beddir

# add CADD scores to SNPs
wait
echo
echo "Adding CADD scores to SNPs..."
bash compileCADDscores.sh ${beddir}
echo "Done compiling CADD scores."
date

echo "Gzipping files..."
parallel --jobs $nproc gzip ::: ${beddir}/individuals/*.bed
parallel --jobs $nproc gzip ::: ${beddir}/individuals/withCADD/*.bed
echo "Done gzipping."
