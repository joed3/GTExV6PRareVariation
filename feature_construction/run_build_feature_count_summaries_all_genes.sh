#!/bin/bash

set -o nounset -o pipefail

dir=${RAREVARDIR}/features
indir_sv=${dir}/variantBeds/individuals

indir=${dir}/bySite
outdir=${dir}/byGene
refdir=${RAREVARDIR}/reference
tss=${refdir}/gencode.v19.genes.v6p.patched_contigs_TSS.bed
gene=${refdir}/gencode.v19.genes.v6p.patched_contigs.bed
exclude=${refdir}/gencode.v19.annotation_coding.lincRNA_padded.bed

# need to run the features before the counts otherwise the features script won't know to build its sub directory structure

# 10 kb
bash build_feature_summaries_all_genes.sh $indir ${outdir}/10kb 10000 $tss
bash build_count_summaries_all_genes.sh $indir $indir_sv ${outdir}/10kb 10000 $tss

# 10 kb without coding variants
bash build_feature_summaries_all_genes.sh $indir ${outdir}/10kb_noPC 10000 $tss $exclude
bash build_count_summaries_all_genes.sh $indir $indir_sv ${outdir}/10kb_noPC 10000 $tss $exclude

# including the gene body
# 10 kb
bash build_feature_summaries_all_genes.sh $indir ${outdir}/10kb_genebody 10000 $gene
bash build_count_summaries_all_genes.sh $indir $indir_sv ${outdir}/10kb_genebody 10000 $gene

# 200 kb
bash build_feature_summaries_all_genes.sh $indir ${outdir}/200kb_genebody 200000 $gene
bash build_count_summaries_all_genes.sh $indir $indir_sv ${outdir}/200kb_genebody 200000 $gene
