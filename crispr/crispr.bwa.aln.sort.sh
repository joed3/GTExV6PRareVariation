#!/usr/bin/bash

# Set number of cores
${ncores} = 10

# Function to align single end read with BWA and sort the output using Samtools. Output is bam file.
alignSort(){
	#fastq=/mnt/lab_data/bassik/gaelenh2/Montgomery_Coding_SNPs/Sequencing_data/GH_SBM_Samp1-66/SBM_Samp${1}_trim.fastq.gz
	fastq=${RAREVARDIR}/data/CRISPR/fastqs/SBM_Samp${1}_trim.fastq.gz
	reference=$2
	mismatches=$3
	outPrefix=${RAREVARDIR}/data/CRISPR/bams/SBM_Samp${1}_n$mismatches

	# Align with BWA
	bwa aln -n $mismatches $reference <(zcat $fastq) > $outPrefix.sai
	bwa samse -f $outPrefix.sam $reference $outPrefix.sai $fastq <(zcat $fastq)

	# Sort and convert to BAM using samtools
	samtools-1.3.1 sort -O BAM -o $outPrefix.bam $outPrefix.sam
	samtools-1.3.1 index $outPrefix.bam $outPrefix.bai

	# Remove intermediate files
	rm $outPrefix.sa*

	# Get summary stats for the alignment results using Picard
	java -jar -Xmx10g ${PICARDPATH}/picard.jar CollectAlignmentSummaryMetrics R=$reference I=$outPrefix.bam O=$outPrefix.summ.stats.txt VALIDATION_STRINGENCY=LENIENT
}

# Define the reference files
cref=<the path to the code subdirectory holding the CRISPR amplicon reference sequences>
gDNAoutlier=$cref/crispr.outlier.gdna.ref.fa
gDNAcontrol=$cref/crispr.control.gdna.ref.fa
cDNAoutlier=$cref/crispr.outlier.cdna.ref.fa
cDNAcontrol=$cref/crispr.control.cdna.ref.fa

# Perform the alignment
export -f alignSort

parallel --jobs ${ncores} alignSort ::: {1..18} ::: $gDNAoutlier ::: 0
parallel --jobs ${ncores} alignSort ::: {19..36} ::: $cDNAoutlier ::: 0
parallel --jobs ${ncores} alignSort ::: {37..51} ::: $gDNAcontrol ::: 0
parallel --jobs ${ncores} alignSort ::: {52..66} ::: $cDNAcontrol ::: 0

