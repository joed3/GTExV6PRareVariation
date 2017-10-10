# GTExV6PRareVariation
Repository to reproduce analyses from the GTEx V6P Rare Variation Manuscript

# To run the code

## Install the dependencies
R packages
* cowplot
* data.table
* doMC
* ggplot2
* matrixStats
* peer
* plotrix
* plyr
* pROC
* RColorBrewer
* reshape2
* scales

Python modules
* numpy
* pybedtools
* pysam
* scipy

Unix packages
* GNU parallel
* tabix

External software
* bedtools v2.26.0 or later
* vcftools
* samtools
* bwa
* PICARD

## Download required files
Download from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz: <br>
* processed data directory (referred to below as \<processed_data\>)

Download from http://www.gtexportal.org/home/datasets: <br>
* `gencode.v19.genes.v6p_model.patched_contigs.gtf.gz`
* `GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz` (gunzip it)
* `GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz` (gunzip it)
* `GTEx_Data_V6_Annotations_SampleAttributesDS.txt`
* `GTEx_Analysis_2015-01-12_Annotations_SubjectPhenotypesDS.txt`
* the covariates used during eQTL discovery
* the eGene and significant variant-gene associations based on permutations for each tissue (`GTEx_Analysis_v6p_eQTL.tar` on the portal)

Download from Gencode (http://www.gencodegenes.org/releases/19.html; Comprehensive gene annotation gtf):
* `gencode.v19.annotation.gtf.gz` (gunzip it)

Download from dbGaP: <br>
(It is possible some of the file names may be different in the final release.)
* `GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz`
* ASE files

Download from http://krishna.gs.washington.edu/download/CADD/v1.2/:
* `whole_genome_SNVs.tsv.gz`
* `whole_genome_SNVs.tsv.gz.tbi`
* `whole_genome_SNVs_inclAnno.tsv.gz`
* `whole_genome_SNVs_inclAnno.tsv.gz.tbi`

Download from ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/:
* `ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`

Download from http://compbio.mit.edu/encode-motifs/:
* `matches.txt.gz`

Download Epigenomics Roadmap files from http://www.broadinstitute.org/~meuleman/reg2map/HoneyBadger2_release/DNase/p2/:
* `prom/BED_files_per_sample/regions_prom_E*.bed`
* `enh/BED_files_per_sample/regions_enh_E*.bed`
* `dyadic/BED_files_per_sample/regions_dyadic_E*.bed`

Download ExAC constraint data from ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/functional_gene_constraint/:
* `forweb_cleaned_exac_r03_march16_z_data_pLI.txt`

Download RNA-Seq data for K562 cell lines from ENCODE at https://www.encodeproject.org/search/?type=Experiment&biosample_term_name=K562&assay_title=RNA-seq
* ENCFF104VTJ_CSHL_1.tsv
* ENCFF201HGA_CSHL_2.tsv
* ENCFF387ZRA_Caltech_1.tsv
* ENCFF870QAL_Caltech_2.tsv
* ENCFF553DDU_UConn_1.tsv
* ENCFF811VBA_UConn_2.tsv
* gencode.v24.tRNAs.gtf.gz

Download the hg19 human reference genome in FASTA format from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz. Combine FASTA files for each chromosome into a single reference FASTA.


Other files: <br>
* `gtex.lumpy.gs.svscore.low_conf.vcf.gz` (Obtained directly from the Hall Lab.)

## Setup
The code relies on an assumed directory structure.
Everything is under an upper level directory.
If you are using the processed files available online, set this upper-level directory to be the path to \<processed_data\>. <br>
Otherwise, set this to be any path where you want scripts to output.

Either run or add the following to your .bashrc (omit trailing slashes for all directories): <br>
```
export RAREVARDIR=<the path to the upper-level directory>
export CADD_DIR=<the path to the CADD files>
export KG_DIR=<the path to the 1000 genomes files>
export ENCODE_MOTIF_DIR=<the path to matches.txt.gz (excluding the file name)>
export ER_DIR=<the path to the Epigenomics Roadmap directory (sub directories are prom, enh, and dyadic)>
export EXAC_DIR=<the path to forweb_cleaned_exac_r03_march16_z_data_pLI.txt (excluding the file name)>
export GTEXCISDIR=<the path to the eGene and significant variant-gene pairs for each tissue from the GTEx v6p release>
export HG19=<the path to the HG19 reference genome in FASTA format>
export PICARDPATH=<the path to the PICARD executable directory>
```

If you are going to run everything from scratch, first make the directories that are assumed to exist. <br>
You can skip this step if you will use the processed data files downloaded from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
```
mkdir ${RAREVARDIR}/logs
mkdir ${RAREVARDIR}/preprocessing
mkdir ${RAREVARDIR}/preprocessing/PEER
mkdir ${RAREVARDIR}/reference
mkdir ${RAREVARDIR}/data
mkdir ${RAREVARDIR}/data/medz
mkdir ${RAREVARDIR}/data/singlez
mkdir ${RAREVARDIR}/data/metasoft
mkdir ${RAREVARDIR}/features
mkdir ${RAREVARDIR}/features/annotations
mkdir ${RAREVARDIR}/features/annotations/ACMG
mkdir ${RAREVARDIR}/features/annotations/ClinVar
mkdir ${RAREVARDIR}/features/annotations/GWAS
mkdir ${RAREVARDIR}/features/annotations/OMIM
mkdir ${RAREVARDIR}/features/annotations/Orphanet
mkdir ${RAREVARDIR}/features/annotations/Other
```

These directories are required to hold the results of scripts for which we do not provide processed data.
```
mkdir ${RAREVARDIR}/paper_figures
mkdir ${RAREVARDIR}/features/variantBeds
mkdir ${RAREVARDIR}/data/CRISPR
mkdir ${RAREVARDIR}/data/CRISPR/bams
mkdir ${RAREVARDIR}/data/K562
```

Then you can run the code below. <br>
Some of the steps reproduce processed files available on https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz and are marked as such.

# Pipeline
## Expression data correction and normalization
Generates processed data that can be downloaded from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>

#### Generate rpkm and read count matrices from GTEx V6P combined file
```
cat <path to downloaded file>/GTEx_Data_V6_Annotations_SampleAttributesDS.txt | \
	cut -f1,14 | sed 's/ - /_/' | sed 's/ /_/g' | sed 's/(//' | sed 's/)//' | sed 's/c-1/c1/' | grep -v NA12878 > \
	${RAREVARDIR}/preprocessing/gtex_2015-01-12_samples_tissues.txt

python preprocessing/split_by_tissues.py \
    --GTEX <path to downloaded file>/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct \
    --SAMPLE ${RAREVARDIR}/preprocessing/gtex_2015-01-12_samples_tissues.txt \
    --OUT ${RAREVARDIR}/preprocessing/PEER \
    --END .rpkm.txt <br>

python preprocessing/split_by_tissues.py \
    --GTEX <path to downloaded file>/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct \
    --SAMPLE ${RAREVARDIR}/preprocessing/gtex_2015-1-12_samples_tissues.txt \
    --OUT ${RAREVARDIR}/preprocessing/PEER \
    --END .reads.txt
```

#### Run bash scripts to generate PEER corrected data (includes non-EAs) with covariates removed
(Uses multiple cores. Currently set to 10 cores. Can be altered by changing the number of cores <br> 
specified by `parallel --jobs 10` in the scripts `preprocessing/PEER/calc.PEER.factors.all.tissues.sh` and <br>
`preprocessing/PEER/calc.residuals.sh`) <br>
```
bash preprocessing/PEER/PEER.pipeline.sh
```

#### Make list of individuals and tissues
```
bash preprocessing/get_tissue_by_individual.sh
```

#### Make flat files from raw RPKMs and PEER-corrected data for all tissues and individuals
```
python preprocessing/gather_filter_normalized_expression.py
python preprocessing/gather_filter_rpkm.py
```

#### Make list of expressed genes
```
cat preprocessing/PEER/*peer.ztrans.txt | cut -f1 | sort | uniq | grep -v Id > preprocessing/gtex.expressed.genes.txt
```

#### Make summary statistics for expressed genes
(Uses multiple cores. Can set the number at the top of the script. Currently set to 10 cores.)
```
Rscript preprocessing/rpkm.expression.analysis.R
```

## Preparing reference files used later
Generates processed data that can be downloaded from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>

Copy or move the GTEx annotation file (`gencode.v19.genes.v6p_model.patched_contigs.gtf.gz`) to `${RAREVARDIR}/reference`.
```
bash preprocessing/process.reference.files.sh <path to>/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz <path to>/gencode.v19.annotation.gtf
```
(relies on `pad.gtf.exons.py`, `gtf2TSS.sh`, and `gtf2genebed.sh`)


## Outlier calling
Generates processed data that can be downloaded from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>

#### Call multi-tissue outliers
(Uses 12 cores. Can set the number at the top of the script.)
```
Rscript call_outliers/call_outliers_medz.R
```

#### Call single-tissue outliers
```
python call_outliers/call_outliers_single_tissue.py
```

#### Compare single-tissue and multi-tissue outliers as well as get stats on each
```
Rscript call_outliers/compare_single_multi_outliers.R
```

#### Run replication for single-tissue and multi-tissue outliers
(The multi-tissue replication uses 12 cores. Can set the number at the top of the script.)
```
Rscript call_outliers/multi_tissue_replication.R
Rscript call_outliers/single_tissue_replication.R
```

## Feature generation

#### Processing VCFs into bed files for each individual
```
bash feature_construction/vcf2bedfiles.sh \
	 <path to>/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz \
	 <path to>/gtex.lumpy.gs.svscore.low_conf.vcf.gz
```
(this script and some of its dependencies use multiple cores [set number at top of relevant scripts]; relies on :
* `vcf2bedfiles_helper_processVCF.sh`
* `vcf2bedfiles_helper_processVCF_SV.sh`
* `vcf2bedfiles_helper_processVCFtoolsOutput.sh`
* `vcf2bedfiles_helper_processVCFtoolsOutput_CNV.sh`
* `compileCADDscores.sh`
* `extractCADDscores_ekt.py`)

#### Extract features to be combined with the individual bed files
```
bash feature_construction/extract.1kg.AF.sh
```
(uses 15 cores, the number of which is set at the top of the script; relies on `process.1kg.AF.py`)
```
bash feature_construction/subset.CADD.features.sh
```
(uses 8 cores, set in the sort command)
```
bash feature_construction/TFBS_pipeline.sh
```
(relies on `pouya.raw.summary.py`)
```
bash feature_construction/ER_pipeline.sh
```

#### Add extracted features to individual bed files
```
bash run_add_features_variant_beds.sh
```
**Important:** Make sure to use bedtools version 2.26.0 or later.
Memory leak in previous versions causes the memory for this script to blow up. <br>
(uses 15+ cores. Set the number of processes at the top of the script. Relies on `add_features_variant_beds.sh`.)

#### Collapse site-level features created above into gene-level features
```
bash feature_construction/run_build_feature_count_summaries_all_genes.sh
```
(uses multiple cores; relies on :
* `build_count_summaries_all_genes.sh` set number of processes at top of script
* `build_feature_summaries_all_genes.sh` set number of processes at top of script
* `build_feature_set.py`)

#### Compile features for outliers and controls
```
bash feature_construction/run_compile_features_outliers.sh
```
(uses up to 10 cores; relies on:
* `compile_features_outliers.sh` set number of processes at top of script
* `compile_features_outliers_nothresh.sh`
* `compile_features_outliers_singletissue.sh` set number of processes at top of script
* `pick_outliers_controls_imbalanced.py`)

## Disease gene annotations
We are providing the processed gene lists for the eight disease gene sets we analyzed for overlap with genes with multi-tissue outliers. <br>
We are also providing, where applicable, the commands and raw files needed to generate these processed lists.

#### ACMG
Source: http://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/ <br>
Raw file: `acmg.csv` <br>
Download from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move the downloaded file to `${RAREVARDIR}/features/annotations/ACMG/`.

#### ClinVar
Source: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id <br>
Raw file: `gene_condition_source_id` <br>
Download from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move the downloaded file to `${RAREVARDIR}/features/annotations/ClinVar/`.

#### GWAS
Source: http://www.ebi.ac.uk/gwas/ <br>
Raw file: `gwas_catalog_v1.0-downloaded_2015-11-30.tsv` <br>
Download from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move the downloaded file to `${RAREVARDIR}/features/annotations/GWAS/`.

#### OMIM
Source: http://www.omim.org/ <br>
Raw files: `morbidmap.txt` and `mim2gene.txt` <br>
Processed file: `omim.genes.txt`
Download the raw and processed files from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move these files to `${RAREVARDIR}/features/annotations/OMIM/`.

To produce the processed file from the raw files: <br>
```
grep '(3)' ${RAREVARDIR}/features/annotations/OMIM/morbidmap.txt | cut -f2 | sed 's/, /\n/g' | sort | uniq > ${RAREVARDIR}/features/annotations/OMIM/omim.genes.temp.txt

grep -wf ${RAREVARDIR}/features/annotations/OMIM/omim.genes.temp.txt ${RAREVARDIR}/features/annotations/OMIM/mim2gene.txt  > ${RAREVARDIR}/features/annotations/OMIM/temp.mim2gene.intersection.txt

cut -f4,5 ${RAREVARDIR}/features/annotations/OMIM/temp.mim2gene.intersection.txt | sort | uniq > ${RAREVARDIR}/features/annotations/OMIM/omim.genes.txt

rm ${RAREVARDIR}/features/annotations/OMIM/omim.genes.temp.txt ${RAREVARDIR}/features/annotations/OMIM/temp.mim2gene.intersection.txt
```

#### Orphanet
Source: http://www.orphadata.org/data/xml/en_product6.xml <br>
Raw file: `en_product6.xml` <br>
Processed file: `orphanet.genes.txt` <br>
Download the raw and processed files from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move these files to `${RAREVARDIR}/features/annotations/Orphanet/`.

To produce the processed file from the raw file: <br>
```
grep ENSG ${RAREVARDIR}/features/annotations/Orphanet/en_product6.xml | sort | uniq | grep -o 'ENSG[0-9]*' > ${RAREVARDIR}/features/annotations/Orphanet/orphanet.genes.txt
```

#### DDG2P
Source: http://www.ebi.ac.uk/gene2phenotype/downloads <br>
Raw file: `DDG2P_2_8_2017.csv.gz` <br>
Download the raw file from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move this file to `${RAREVARDIR}/features/annotations/DDG2P/`.

#### Other: Cardiovascular and Cancer disease genes
We assessed overlap of genes with multi-tissue outliers with two expert curated disease gene lists: one for heritable cancer predisposition and one for heritable cardiovascular disease. See the methods section of our manuscript for more information. <br>
Raw files: `cancer.genes.gold.standard.csv` (Cancer), `cardio.genes.gold.standard.csv` (Cardio) 
Download the raw files from https://s3-us-west-2.amazonaws.com/gtex-v6p-rare-variation-data/GTExV6PRareVariationData.tar.gz. <br>
Move these files to `${RAREVARDIR}/features/annotations/Other/`.

## Shared eQTLs defined by METASOFT
Process the METASOFT results choosing the single best variant tested per gene as determined by P-value from the RE2 model. <br>
Also provide summary statistics regarding the number of tissues the eQTL is active in and the expression level for the gene <br>
across tissues.
```
python shared.eqtls/bf.metasoft.py --META ${RAREVARDIR}/data/metasoft/Metasoft_Output_v6p.txt \
    --TISS ${RAREVARDIR}/data/metasoft/Metasoft_tissue_order.txt \
    --OUT ${RAREVARDIR}/data/metasoft/gtex.metasoft.v6p.selected.txt
Rscript shared.eqtls/metasoft.summary.R
```

## Validation of large-effect rare variants using CRISPR-Cas9 genome editing
Prioritize variants for validation with CRISPR
```
Rscript crispr/prioritize.for.crispr.R
```

Prune the list of prioritized CRISPR variants down to a manageable size (N ~ 12)
```
Rscript crispr/prune.crispr.variants.R
```

Add major/minor allele info
```
python crispr/add.major.minor.alleles.py --IN ${RAREVARDIR}/data/CRISPR/crispr.overexpression.candidates.pruned.vcf \
    --OUT ${RAREVARDIR}/data/CRISPR/crispr.overexpression.candidates.pruned.maf.alleles.vcf \
    --FRQ ${RAREVARDIR}/features/variantBeds/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller_123EAonly_SNPs.frq
```

Extract 100 bp sequences from the hg19 reference centered on each variant position
```
out=${RAREVARDIR}/data/CRISPR/crispr.candidates.donor.seq.raw.fa
if [ -e $out ]; then
    rm $out
fi
ref=${HG19}/hg19.fa
positions=`tail -n +2 ${RAREVARDIR}/data/CRISPR/crispr.candidates.pruned.vep.loftee.parsed.vcf | awk '{print "chr"$1":"$2-49"-"$2+49}'`
for line in $positions; do
    samtools faidx $ref $line >> $out
done
```

Generate donor sequences. One sequence for wild-type and one for rare allele.
```
python crispr/process.crispr.donor.seq.py \
    --VCF ${RAREVARDIR}/data/CRISPR/crispr.candidates.pruned.vep.loftee.parsed.vcf \
    --FASTA ${RAREVARDIR}/data/CRISPR/crispr.candidates.donor.seq.raw.fa \
    --OUT ${RAREVARDIR}/data/CRISPR/crispr.candidates.donor.seq.fa
```

Index the reference FASTA files for each amplicon region in the cDNA and gDNA
```
bwa index -p crispr/crispr.outlier.cdna.ref.fa crispr/crispr.outlier.cdna.ref.fa
bwa index -p crispr/crispr.control.cdna.ref.fa crispr/crispr.control.cdna.ref.fa
bwa index -p crispr/crispr.outlier.gdna.ref.fa crispr/crispr.outlier.gdna.ref.fa
bwa index -p crispr/crispr.control.gdna.ref.fa crispr/crispr.control.gdna.ref.fa
```

Perform mapping with BWA, sort and output to BAM format using samtools. <br>
Uses 10 cores. This is set at the top of the script.
```
bash crispr/crispr.bwa.aln.sort.sh
```

Summarize the CRISPR results
```
R CMD BATCH --no-save crispr/summarize.crispr.results.R
```

## Main figures
#### Figure 1
```
bash paper_figures/pick.cartoon.example.sh
Rscript paper_figures/figure1a.plot.cartoon.example.R
Rscript paper_figures/figure1b.outlier.sharing.R
Rscript paper_figures/figure1c.replication.rate.consistent.R
```

#### Figure 2
```
bash paper_figures/count_rarevars.sh
Rscript paper_figures/figure2a.count.enrichments.R
Rscript paper_figures/figure2b.threshold.enrichments.R
Rscript paper_figures/figure2c.ASE.enrichments.R
Rscript paper_figures/figure2.R
```

#### Figure 3
```
Rscript paper_figures/figure3a.rare.variant.class.enrichments.R
Rscript paper_figures/figure3b.feature.enrichments.R
Rscript paper_figures/figure3de.outlier.effect.size.R
Rscript paper_figures/figure3.R
```

#### Figure 4
```
Rscript paper_figures/figure4a.uk10k.R
Rscript paper_figures/figure4b.exac.enrichments.R
Rscript paper_figures/figure4c.gene.list.enrichments.R
Rscript paper_figures/figure4.R
```

#### Figure 5
```
Rscript getPosteriorsEval.R data/genomic_features.txt data/outliers.txt
Rscript paper_figures/figure5b.R
Rscript getPosteriorsApp.R data/genomic_features.txt data/outliers.txt data/postprobs_all.txt
Rscript paper_figures/figure5c.R
```

## Supplemental figures

#### GTEx design matrix and the MAFs in 1KG european populations
(Running `calculate.euro.subpop.af.sh` uses five cores. This is set at the top of the script.)
```
bash feature_construction/calculate.euro.subpop.af.sh
Rscript feature_construction/assess.euro.subpop.R
Rscript paper_figures/suppfig.gtex.design.R
```

#### Number of rare variants per individual and PCA
You need to set the path to the downloaded expression covariates and subject annotations in the script below. <br> 
```
Rscript paper_figures/suppfig.number.rare.vars.pca.R
```

#### Improvement of replication of outliers across tissues by PEER correction
```
Rscript paper_figures/ExtendedDataFigure3.R
```

#### Distribution of the number of genes with a multi-tissue outlier
You need to set the path to the subject annotations in the script below. <br> 
```
Rscript paper_figures/suppfig.number.outliers.per.individual.R
```

#### Single-tissue replication analysis controlling for sampling differences
Generated when running `paper_figures/figure1b.outlier.sharing.R` above

#### Overlap between single and multi-tissue outliers
```
Rscript paper_figures/suppfig.compare.single.multi.R
```

#### Single-tissue outlier rare variant enrichments
Generated when running `paper_figures/figure2b.threshold.enrichments.R` above

#### Multi-tissue outlier rare variant enrichments when excluding coding regions
Generated when running `paper_figures/figure2a.count.enrichments.R above`

#### Comparison of enrichment of functional rare SNVs between over- and underexpression outliers
Generated when running `paper_figures/figure3b.feature.enrichments.R`

#### Enrichment of functional genomic annotations among an expanded set of multi-tissue outliers
```
Rscript paper_figures/EDF7.R
```

#### Comparing distribution of the number of rare variants in the GTEx cohort for disease and control genes
```
bash paper_figures/annotate.variants.by.gene.sh 
Rscript paper_figures/suppfig.rare.var.counts.disease.genes.gtex.cohort.R
```

#### Evolutionary constraint of genes with multi-tissue outliers and their overlap with multi-tissue eGenes
Relies on the eGene and singificant variant-gene associations as downloaded for the GTEx v6p release from the portal. Also relies on output from `paper_figures/figure4b.exac.enrichments.R`, `paper_figures/figure4c.gene.list.enrichments.R`, and `paper_figures/suppfig.rare.var.counts.disease.genes.gtex.cohort.R`.
```
Rscript paper_figures/suppfig.egene.enrichment.R
```

#### Association between ASE and RIVER scores
```
Rscript paper_figures/EDF9d.R
```

#### Correlation between test posterior probabilities and the fraction of tissues
```
Rscript paper_figures/main_RIVER_10CV.R
Rscript paper_figures/EDF9c.R
```

#### Distribution of predictive scores for pathogenic variants and all variants
```
Rscript paper_figures/EDF10bd.R
```

#### Expression levels for genes proximal to pathogenic variants
```
Rscript paper_figures/EDF10ef.rpkm.R
Rscript paper_figures/EDF10ef.Zscores.R
```
#### Adjusted R-squared values between top 15 PEER factors and top 20 sample and subject covariates in skeletal muscle
``` 
Rscript paper_figures/muscle_covariates_peerfactors.R
Rscript paper_figures/superheat_peer_muscle.R /data/muscle_samples_covariates_peerFactors.RData muscle.sample.peer.pdf 0.7
Rscript paper_figures/superheat_peer_muscle.R /data/muscle_subject_covariates_peerFactors.RData muscle.subject.peer.pdf 0.7
```

#### Adjusted R-squared values between the total expression component removed by PEER in each of 44 tissues and top 20 sample and subject covariates
``` 
Rscript paper_figures/pve_samples_pertiss.R
Rscript paper_figures/pve_subject_pertiss.R
Rscript paper_figures/process_results.R
Rscript paper_figures/superheat_expression.R /data/superheat.subject.RData rv.subject.expression.pdf 0.25
Rscript paper_figures/superheat_expression.R /data/superheat.sample.RData rv.sample.expression.pdf 0.25
```

#### Comparison between the predictive power of RIVER and that of the genomic annotation model, as in Fig. 5a, across different Z-score thresholds for outlier calling
``` 
Rscript paper_figures/main_RIVER_VaryingThrds.R
Rscript paper_figures/Generate_figures_RIVER_VaryingThrds.R
```
