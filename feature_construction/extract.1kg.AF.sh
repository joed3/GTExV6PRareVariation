#!/bin/bash

set -o nounset

# get AF of the 5 superpopulations from the 1kg vcfs

# for each chromosome, extract the AF information from the vcf
# only do this for the sites in GTEx

dir=${RAREVARDIR}/features/variantBeds
gtexsnp=${dir}/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller_123EAonly_SNPs.frq
gtexindel=${dir}/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller_123EAonly_indels.frq

vcfdir=$KG_DIR
outdir=${dir}/1KG

if [ ! -d $outdir ]; then
    mkdir $outdir
fi

export gtexsnp
export gtexindel
export vcfdir
export outdir

for i in chr{1..22}; do echo $i; done | \
parallel --jobs 15 "vcftools --gzvcf ${vcfdir}/ALL.{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --out ${outdir}/{}.indels --keep-only-indels --positions $gtexindel --get-INFO EAS_AF --get-INFO AMR_AF --get-INFO AFR_AF --get-INFO EUR_AF --get-INFO SAS_AF"

for i in chr{1..22}; do echo $i; done | \
parallel --jobs 15 "vcftools --gzvcf ${vcfdir}/ALL.{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --out ${outdir}/{}.SNPs --remove-indels --positions $gtexsnp --get-INFO EAS_AF --get-INFO AMR_AF --get-INFO AFR_AF --get-INFO EUR_AF --get-INFO SAS_AF"

cat ${outdir}/*indels.INFO | ./process.1kg.AF.py | sort -k1,1 -k2,2n > ${outdir}/indels.1kg.AF.bed

cat ${outdir}/*SNPs.INFO | ./process.1kg.AF.py | sort -k1,1 -k2,2n > ${outdir}/SNPs.1kg.AF.bed

rm ${outdir}/chr*.INFO
gzip ${outdir}/*.bed
