#!/bin/bash

set -o nounset -o errexit -o pipefail

ncores=5

# calculate the 1KG AF for each of the european subpopulations

# sample info file downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/
# and provided in this repository

outdir=${RAREVARDIR}/features/variantBeds/1KG
prefix="${KG_DIR}/ALL.chr"
suffix=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
samplefile=20130606_sample_info.txt
vcfinds=${outdir}/individuals.in.vcf.txt

# extract individuals in phase 3 from the vcf
zcat ${prefix}10${suffix} | head -n 5000 | grep CHROM | cut --complement -f1-9 | tr '\t' '\n' > $vcfinds

# function to actually run the AF calculation
subpop_af() {
subpop=$1
keepfile=${outdir}/${subpop}.inds.txt

# get individuals to keep for the given subpopulation 
# (first subsetting to individuals in the vcf)
cat $vcfinds | awk '{print "^"$1}' | grep -f - $samplefile | \
awk -v subpop=$subpop 'BEGIN{FS="\t"}{if($3==subpop){print $1}}' > $keepfile

# chromosome 1 (to overwrite output files if the exist, and have the header)
vcftools --gzvcf ${prefix}1${suffix} --stdout --remove-indels --keep $keepfile --freq2 > ${outdir}/${subpop}.SNPs.AF.txt
vcftools --gzvcf ${prefix}1${suffix} --stdout --keep-only-indels --keep $keepfile --freq2 > ${outdir}/${subpop}.indels.AF.txt

# all others
for chrom in {2..22}
do
    vcf="$prefix$chrom$suffix"
    vcftools --gzvcf $vcf --stdout --remove-indels --keep $keepfile --freq2 | tail -n +2 >> ${outdir}/${subpop}.SNPs.AF.txt
    vcftools --gzvcf $vcf --stdout --keep-only-indels --keep $keepfile --freq2 | tail -n +2 >> ${outdir}/${subpop}.indels.AF.txt
done
}

# actually run things
export -f subpop_af
export outdir
export prefix
export suffix
export samplefile
export vcfinds
parallel --jobs $ncores subpop_af ::: CEU TSI FIN GBR IBS
