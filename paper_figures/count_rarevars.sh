#!/bin/bash

set -o nounset -o errexit -o pipefail

# count the number of non-coding rare variants genome-wide as well as all variants

# non-coding here is defined as being in outside of lincRNA and protein-coding gene exons (internal exons padded by 5 bp)
# note that this also excludes UTRs
dir=${RAREVARDIR}
coding=${dir}/reference/gencode.v19.annotation_coding.lincRNA_padded.bed
ind_file=${dir}/preprocessing/gtex_2015-01-12_wgs_ids.txt #121
snv_indel_dir=${dir}/features/bySite
sv_dir=${dir}/features/variantBeds/individuals
outfile=${dir}/data/variant_counts_per_individual_all.txt
outfile2=${dir}/data/variant_counts_per_individual_noncoding.txt
outfilekg=${dir}/data/variant_counts_per_individual_filtered_1KG.txt

echo -e "Ind\tSNV\tindel\tSV" > $outfile
echo -e "Ind\tSNV\tindel\tSV" > $outfile2
echo -e "Ind\tSNV\tindel" > $outfilekg

# all
for sample in `cat $ind_file`
do
    snv_count=`zcat ${snv_indel_dir}/${sample}_SNPs_features.bed.gz | awk '$4<=0.01 && ($(NF-1)=="NA" || $(NF-1)<=0.01) {print $1"\t"$2"\t"$3}' | wc -l`
    indel_count=`zcat ${snv_indel_dir}/${sample}_indels_features.bed.gz | awk '$4<=0.01 && ($(NF-1)=="NA" || $(NF-1)<=0.01) {print $1"\t"$2"\t"$3}' | wc -l`
    # we don't have SV calls for one individual
    if [ $sample == "GTEX-WHWD" ]; then
	sv_count="NA"
    else
	sv_count=`zcat ${sv_dir}/${sample}_HallLabSV.bed.gz | awk '$4<=0.01' | wc -l`
    fi
    echo -e "$sample\t$snv_count\t$indel_count\t$sv_count" >> $outfile
done

# noncoding
for sample in `cat $ind_file`
do
    snv_count=`zcat ${snv_indel_dir}/${sample}_SNPs_features.bed.gz | awk '$4<=0.01 && ($(NF-1)=="NA" || $(NF-1)<=0.01) {print $1"\t"$2"\t"$3}' | bedtools intersect -sorted -v -a stdin -b $coding | wc -l`
    indel_count=`zcat ${snv_indel_dir}/${sample}_indels_features.bed.gz | awk '$4<=0.01 && ($(NF-1)=="NA" || $(NF-1)<=0.01) {print $1"\t"$2"\t"$3}' | bedtools intersect -sorted -v -a stdin -b $coding | wc -l`
    # we don't have SV calls for one individual
    if [ $sample == "GTEX-WHWD" ]; then
	sv_count="NA"
    else
	sv_count=`zcat ${sv_dir}/${sample}_HallLabSV.bed.gz | awk '$4<=0.01' | bedtools intersect -sorted -v -a stdin -b $coding | wc -l`
    fi
    echo -e "$sample\t$snv_count\t$indel_count\t$sv_count" >> $outfile2
done

# count the number of rare variants exlcuded by the 1KG filter (not relevant for SVs)
for sample in `cat $ind_file`
do
    snv_count=`zcat ${snv_indel_dir}/${sample}_SNPs_features.bed.gz | awk '$4<=0.01 && $(NF-1)!="NA" && $(NF-1)>0.01' | wc -l`
    indel_count=`zcat ${snv_indel_dir}/${sample}_indels_features.bed.gz | awk '$4<=0.01 && $(NF-1)!="NA" && $(NF-1)>0.01' | wc -l`
    echo -e "$sample\t$snv_count\t$indel_count" >> $outfilekg
done


