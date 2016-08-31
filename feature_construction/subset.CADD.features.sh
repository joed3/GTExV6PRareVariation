#!/bin/bash

# Chrom, Pos0, Pos1, CpG, priPhCons, mamPhCons, verPhCons, priPhyloP, mamPhyloP, verPhyloP, GerpN, GerpS, GerpRS, GerpRSpval, fitCons, TFBS, TFBSPeaks, TFBSPeaksMax

zcat $CADD_DIR/whole_genome_SNVs_inclAnno.tsv.gz | tail -n+3 | \
	awk 'BEGIN{OFS="\t"}{print "chr"$1,$2-1,$2,$15,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$39,$79,$80,$81}' | uniq | \
	sort -k1,1 -k2,2n --parallel=8 | bgzip -c > ${CADD_DIR}/sorted.consolidated.annotations.bed.gz

tabix -p bed ${CADD_DIR}/sorted.consolidated.annotations.bed.gz
