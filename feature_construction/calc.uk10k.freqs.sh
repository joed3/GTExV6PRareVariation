#!/bin/bash

set -o nounset -o errexit -o pipefail

mkdir ${RAREVARDIR}/features/variantBeds/UK10K
outdir=${RAREVARDIR}/features/variantBeds/UK10K

vcftools --gzvcf <path to>/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz \
	--freq --stdout > ${outdir}/GTEx.v6p.WGS148.variant.sites.txt 

vcf-concat `ls -f <path to>/*.TWINSUK.beagle.anno.csq.shapeit.20131101.vcf.gz` | \
	vcftools --vcf - --positions ${outdir}/GTEx.v6p.WGS148.variant.sites.txt \
	--recode --stdout | \
	bgzip > ${outdir}/GTEx.v6p.WGS148.variant.sites.UK10K.TWINSUK.vcf.gz

vcf-concat `ls -f <path to>/*.ALSPAC.beagle.anno.csq.shapeit.20131101.vcf.gz` | \
	vcftools --vcf - --positions ${outdir}/GTEx.v6p.WGS148.variant.sites.txt \
	--recode --stdout | \
	bgzip > ${outdir}/GTEx.v6p.WGS148.variant.sites.UK10K.ALSPAC.vcf.gz

bcftools merge ${outdir}/GTEx.v6p.WGS148.variant.sites.UK10K.ALSPAC.vcf.gz \
	${outdir}/GTEx.v6p.WGS148.variant.sites.UK10K.TWINSUK.vcf.gz | \
	bgzip > ${outdir}/GTEx.v6p.WGS148.variant.sites.UK10K.TWINSUK.ALSPAC.merged.vcf.gz

vcftools --gzvcf ${outdir}/GTEx.v6p.WGS148.variant.sites.UK10K.TWINSUK.ALSPAC.merged.vcf.gz \
	--freq --keep $1 --stdout > ${outdir}/gtex_allrare.txt
