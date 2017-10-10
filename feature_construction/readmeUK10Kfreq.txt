

# /home/xli6/tools/vcftools_0.1.12b/bin/vcftools --gzvcf /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/genotypes/WGS/WGS148/variant_calls/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz --freq --stdout > /srv/gs1/projects/montgomery/xin/gtex/GENOTYPE/processed/VCFtools/list/GTEx.WGS148.epigenome.GWAS.list.11042015.txt 

# /srv/gs1/projects/montgomery/xin/tools/vcftools_0.1.12b/bin/vcf-concat `ls -f /srv/gs1/projects/montgomery/xin/uk10k/Whole_genome_cohorts_4000/UK10K_COHORT_TWINSUK/VCF/EGAD00001000741/*.TWINSUK.beagle.anno.csq.shapeit.20131101.vcf.gz` | /srv/gs1/projects/montgomery/xin/tools/vcftools_0.1.12b/bin/vcftools --vcf - --positions /srv/gs1/projects/montgomery/xin/gtex/GENOTYPE/processed/VCFtools/list/GTEx.WGS148.epigenome.GWAS.list.11042015.txt --recode --stdout | /srv/gs1/projects/montgomery/xin/tools/tabix-0.2.6/bgzip > /srv/gs1/projects/montgomery/xin/uk10k/VCFnotes/XinProcess/GTExSites/chrAll.TWINSUK.beagle.anno.csq.shapeit.20131101.vcf.gz.GTEx.WGS148.epigenome.GWAS.list.11042015.txt.vcf.gz


# /srv/gs1/projects/montgomery/xin/tools/vcftools_0.1.12b/bin/vcf-concat `ls -f /srv/gs1/projects/montgomery/xin/uk10k/Whole_genome_cohorts_4000/UK10K_COHORT_ALSPAC/VCF/EGAD00001000740/*.ALSPAC.beagle.anno.csq.shapeit.20131101.vcf.gz` | /srv/gs1/projects/montgomery/xin/tools/vcftools_0.1.12b/bin/vcftools --vcf - --positions /srv/gs1/projects/montgomery/xin/gtex/GENOTYPE/processed/VCFtools/list/GTEx.WGS148.epigenome.GWAS.list.11042015.txt --recode --stdout | /srv/gs1/projects/montgomery/xin/tools/tabix-0.2.6/bgzip > /srv/gs1/projects/montgomery/xin/uk10k/VCFnotes/XinProcess/GTExSites/chrAll.ALSPAC.beagle.anno.csq.shapeit.20131101.vcf.gz.GTEx.WGS148.epigenome.GWAS.list.11042015.txt.vcf.gz


# /srv/gs1/projects/montgomery/xin/tools/vcftools_0.1.12b/bin/vcf-merge chrAll.ALSPAC.beagle.anno.csq.shapeit.20131101.vcf.gz.GTEx.WGS148.epigenome.GWAS.list.11042015.txt.vcf.gz chrAll.TWINSUK.beagle.anno.csq.shapeit.20131101.vcf.gz.GTEx.WGS148.epigenome.GWAS.list.11042015.txt.vcf.gz | /srv/gs1/projects/montgomery/xin/tools/tabix-0.2.6/bgzip > /srv/gs1/projects/montgomery/xin/uk10k/VCFnotes/XinProcess/GTExSites/chrAll.TWINSUK.ALSPAC.GTEx.WGS148.epigenome.GWAS.list.11042015.txt.vcf.gz


# /srv/gs1/projects/montgomery/xin/tools/bcftools-1.2/bcftools merge chrAll.ALSPAC.beagle.anno.csq.shapeit.20131101.vcf.gz.GTEx.WGS148.epigenome.GWAS.list.11042015.txt.vcf.gz chrAll.TWINSUK.beagle.anno.csq.shapeit.20131101.vcf.gz.GTEx.WGS148.epigenome.GWAS.list.11042015.txt.vcf.gz | /srv/gs1/projects/montgomery/xin/tools/tabix-0.2.6/bgzip > /srv/gs1/projects/montgomery/xin/uk10k/VCFnotes/XinProcess/GTExSites/chrAll.TWINSUK.ALSPAC.GTEx.WGS148.epigenome.GWAS.list.11042015.bcftools.vcf.gz


# /home/xli6/tools/vcftools_0.1.12b/bin/vcftools --gzvcf /srv/gs1/projects/montgomery/xin/uk10k/VCFnotes/XinProcess/GTExSites/chrAll.TWINSUK.ALSPAC.GTEx.WGS148.epigenome.GWAS.list.11042015.bcftools.vcf.gz --freq --keep $1 --stdout > gtex_allrare.txt



