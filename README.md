# GTExV6PRareVariation
Repository to reproduce analyses from the GTEx V6P Rare Variation Manuscript

# To run the code
## Download required files
Download from \<website - coming soon\>: <br>
* processed data directory (referred to below as \<processed_data\>)

Download from http://www.gtexportal.org/home/datasets: <br>
* gencode.v19.genes.v6p_model.patched_contigs.gtf.gz
* GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz (gunzip it) <br>
* GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz (gunzip it) <br>
* GTEx_Data_V6_Annotations_SampleAttributesDS.txt
* GTEx_Analysis_2015-01-12_Annotations_SubjectPhenotypesDS.txt
* the covariates used during eQTL discovery

Download from Gencode (http://www.gencodegenes.org/releases/19.html; Comprehensive gene annotation gtf):
* gencode.v19.annotation.gtf.gz (gunzip it)

Download from dbGaP: <br>
(It is possible some of the file names may be different in the final release.)
* GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz
* ASE files

Download from http://krishna.gs.washington.edu/download/CADD/v1.2/:
* whole_genome_SNVs.tsv.gz
* whole_genome_SNVs.tsv.gz.tbi
* whole_genome_SNVs_inclAnno.tsv.gz
* whole_genome_SNVs_inclAnno.tsv.gz.tbi

Download from ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/:
* ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

Download from http://compbio.mit.edu/encode-motifs/:
* matches.txt.gz

Download Epigenomics Roadmap files from http://www.broadinstitute.org/~meuleman/reg2map/HoneyBadger2_release/DNase/p2/:
* prom/BED_files_per_sample/regions_prom_E*.bed
* enh/BED_files_per_sample/regions_enh_E*.bed
* dyadic/BED_files_per_sample/regions_dyadic_E*.bed

Download ExAC constraint data from ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/functional_gene_constraint/:
* forweb_cleaned_exac_r03_march16_z_data_pLI.txt

Other files: <br>
* gtex.lumpy.gs.svscore.low_conf.vcf.gz (Obtained directly from the Hall Lab.)

## Setup
The code relies on an assumed directory structure.
Everything is under an upper level directory.
If you are using the processed files available online, set this upper-level directory to be the path to \<processed_data\>. <br>
Otherwise, set this to be any path where you want scripts to output.

Either run or add the following to your .bashrc (omit trailing slashes for all directories): <br>
export RAREVARDIR=\<the path to the upper-level directory\> <br>
export CADD_DIR=\<the path to the CADD files\> <br>
export KG_DIR=\<the path to the 1000 genomes files\> <br>
export ENCODE_MOTIF_DIR=\<the path to matches.txt.gz (excluding the file name)\> <br>
export ER_DIR=\<the path to the Epigenomics Roadmap directory (sub directories are prom, enh, and dyadic)\> <br>
export EXAC_DIR=\<the path to forweb_cleaned_exac_r03_march16_z_data_pLI.txt (excluding the file name)\>

If you are going to run everything from scratch, first make the directories that are assumed to exist. <br>
You can skip this step if you will use the processed data files downloaded from \<website - coming soon\>. <br>
mkdir ${RAREVARDIR}/logs <br>
mkdir ${RAREVARDIR}/preprocessing <br>
mkdir ${RAREVARDIR}/preprocessing/PEER <br>
mkdir ${RAREVARDIR}/reference <br>
mkdir ${RAREVARDIR}/data <br>
mkdir ${RAREVARDIR}/data/medz <br>
mkdir ${RAREVARDIR}/data/singlez <br>

These directories are required to hold the results of scripts for which we do not provide processed data.
mkdir ${RAREVARDIR}/paper_figures <br>
mkdir ${RAREVARDIR}/features <br>
mkdir ${RAREVARDIR}/features/variantBeds <br>
mkdir ${RAREVARDIR}/features/annotations <br>

Then you can run the code below. <br>
Some of the steps reproduce processed files available on \<website - coming soon\> and are marked as such.

# Pipeline
## Expression data correction and normalization
Generates processed data that can be downloaded from \<website - coming soon\>. <br>

#### Generate rpkm and read count matrices from GTEx V6P combined file
cat \<path to downloaded file\>/GTEx_Data_V6_Annotations_SampleAttributesDS.txt | \ <br>
	cut -f1,14 | sed 's/ - /_/' | sed 's/ /_/g' | sed 's/(//' | sed 's/)//' | sed 's/c-1/c1/' | grep -v NA12878 > \ <br>
	${RAREVARDIR}/preprocessing/gtex_2015-01-12_samples_tissues.txt

python preprocessing/split_by_tissues.py \ <br>
    --GTEX \<path to downloaded file\>/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct \ <br>
    --SAMPLE ${RAREVARDIR}/preprocessing/gtex_2015-01-12_samples_tissues.txt \ <br>
    --OUT ${RAREVARDIR}/preprocessing/PEER \ <br>
    --END .rpkm.txt <br>

python preprocessing/split_by_tissues.py \ <br>
    --GTEX \<path to downloaded file\>/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct \ <br>
    --SAMPLE ${RAREVARDIR}/preprocessing/gtex_2015-1-12_samples_tissues.txt \ <br>
    --OUT ${RAREVARDIR}/preprocessing/PEER \ <br>
    --END .reads.txt

#### Run bash scripts to generate PEER corrected data (includes non-EAs) with covariates removed
(Uses multiple cores. Currently set to 10 cores. Can be altered by changing the number of cores <br> 
specified by `parallel --jobs 10` in the scripts `preprocessing/PEER/calc.PEER.factors.all.tissues.sh` and <br>
`preprocessing/PEER/calc.residuals.sh`) <br>
bash preprocessing/PEER/PEER.pipeline.sh

#### Make list of individuals and tissues
bash preprocessing/get_tissue_by_individual.sh

#### Make flat files from raw RPKMs and PEER-corrected data for all tissues and individuals
python preprocessing/gather_filter_normalized_expression.py <br>
python preprocessing/gather_filter_rpkm.py

#### Make list of expressed genes
cat preprocessing/PEER/*peer.ztrans.txt | cut -f1 | sort | uniq | \ <br>
grep -v Id > preprocessing/gtex.expressed.genes.txt

#### Make summary statistics for expressed genes
(Uses multiple cores. Can set the number at the top of the script. Currently set to 10 cores.) <br>
Rscript preprocessing/rpkm.expression.analysis.R >&logs/rpkm.expression.analysis.Rout

## Preparing reference files used later
Generates processed data that can be downloaded from \<website - coming soon\>. <br>

Copy or move the GTEx annotation file (gencode.v19.genes.v6p_model.patched_contigs.gtf.gz) to ${RAREVARDIR}/reference.

bash preprocessing/process.reference.files.sh \<path to\>/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz \<path to\>/gencode.v19.annotation.gtf <br>
(relies on pad.gtf.exons.py, gtf2TSS.sh, and gtf2genebed.sh)


## Outlier calling
Generates processed data that can be downloaded from \<website - coming soon\>. <br>

#### Call multi-tissue outliers
(Uses multiple cores. Can set the number at the top of the script.) <br>
Rscript call_outliers/call_outliers_medz.R

#### Call single-tissue outliers
python call_outliers/call_outliers_single_tissue.py

#### Compare single-tissue and multi-tissue outliers as well as get stats on each
Rscript call_outliers/compare_single_multi_outliers.R

#### Run replication for single-tissue and multi-tissue outliers
(The multi-tissue replication uses multiple cores.)
Rscript call_outliers/multi_tissue_replication.R
Rscript call_outliers/single_tissue_replication.R


## Feature generation

#### Processing VCFs into bed files for each individual
bash feature_construction/vcf2bedfiles.sh 
	 \<path to \>/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz 
	 \<path to\>/gtex.lumpy.gs.svscore.low_conf.vcf.gz <br>
(relies on :
* vcf2bedfiles_helper_processVCF.sh
* vcf2bedfiles_helper_processVCF_SV.sh
* vcf2bedfiles_helper_processVCFtoolsOutput.sh
* vcf2bedfiles_helper_processVCFtoolsOutput_CNV.sh
* compileCADDscores.sh
* extractCADDscores_ekt.py)

#### Extract features to be combined with the individual bed files
bash feature_construction/extract.1kg.AF.sh <br>
(uses multiple cores; relies on process.1kg.AF.py)

bash feature_construction/subset.CADD.features.sh <br>

bash feature_construction/TFBS_pipeline.sh <br>
(relies on pouya.raw.summary.py)

bash feature_construction/ER_pipeline.sh

#### Add extracted features to individual bed files
bash run_add_features_variant_beds.sh <br>
**Important:** Make sure to use bedtools version 2.26.0 or later.
Memory leak in previous versions causes the memory for this script to blow up. <br>
(uses multiple cores; relies on add_features_variant_beds.sh)

#### Collapse site-level features created above into gene-level features
bash feature_construction/run_build_feature_count_summaries_all_genes.sh
(uses multiple cores; relies on :
* build_count_summaries_all_genes.sh
* build_feature_summaries_all_genes.sh
* build_feature_set.py)

#### Compile features for outliers and controls
bash feature_construction/run_compile_features_outliers.sh
(uses multiple cores; relies on:
* compile_features_outliers.sh
* compile_features_outliers_nothresh.sh
* compile_features_outliers_singletissue.sh
* pick_outliers_controls_imbalanced.py)

## Disease gene annotations
We are providing the processed gene lists for the eight disease gene sets we analyzed for overlap with genes with multi-tissue outliers. <br>
We are also providing, where applicable, the commands and raw files needed to generate these processed lists.

#### ACMG
Source: http://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/ <br>
Raw file: disease.genes/ACMG/acmg.csv

#### ClinVar
Source: http://www.ncbi.nlm.nih.gov/clinvar/ <br>
Raw file: disease.genes/ClinVar/gene_condition_source_id

#### GWAS
Source: http://www.ebi.ac.uk/gwas/ <br>
Raw file: disease.genes/GWAS/gwas_catalog_v1.0-downloaded_2015-11-30.tsv

#### OMIM
Source: http://www.omim.org/ <br>
Raw files: disease.genes/OMIM/morbidmap.txt and disease.genes/OMIM/mim2gene.txt <br>
Processed file: disease.genes/OMIM/omim.genes.txt

To produce the processed file from the raw files: <br>
grep '(3)' disease.genes/OMIM/morbidmap.txt | cut -f2 | sed 's/, /\n/g' | sort | uniq > disease.genes/OMIM/omim.genes.temp.txt

grep -wf disease.genes/OMIM/omim.genes.temp.txt disease.genes/OMIM/mim2gene.txt  > disease.genes/OMIM/temp.mim2gene.intersection.txt

cut -f4,5 disease.genes/OMIM/temp.mim2gene.intersection.txt | sort | uniq > disease.genes/OMIM/omim.genes.txt

rm disease.genes/OMIM/omim.genes.temp.txt disease.genes/OMIM/temp.mim2gene.intersection.txt

#### Orphanet
Source: http://www.orpha.net/ <br>
Raw file: disease.genes/Orphanet/en_product6.xml <br>
Processed file: disease.genes/Orphanet/orphanet.genes.txt

To produce the processed file from the raw file: <br>
grep ENSG disease.genes/Orphanet/en_product6.xml | sort | uniq | grep -o 'ENSG[0-9]*' > disease.genes/Orphanet/orphanet.genes.txt

#### Other: Cardiovascular and Cancer disease genes
We assessed overlap of genes with multi-tissue outliers with two expert curated disease gene lists: one for heritable cancer predisposition and one for heritable cardiovascular disease. See the methods section of our manuscript for more information. <br>
Raw files: cancer.genes.gold.standard.csv (Cancer), cardio.genes.gold.standard.csv (Cardio) 

## Main figures
#### Figure 1
bash paper_figures/pick.cartoon.example.sh <br>
Rscript paper_figures/figure1a.plot.cartoon.example.R <br>
Rscript paper_figures/figure1b.outlier.sharing.R <br>
Rscript paper_figures/figure1c.replication.rate.consistent.R

#### Figure 2
bash paper_figures/count_rarevars.sh <br>
Rscript paper_figures/figure2a.count.enrichments.R <br>
Rscript paper_figures/figure2b.threshold.enrichments.R <br>
Rscript paper_figures/figure2c.ASE.enrichments.R <br>
Rscript paper_figures/figure2.R

#### Figure 3
Rscript paper_figures/figure3a.rare.variant.class.enrichments.R <br>
Rscript paper_figures/figure3b.feature.enrichments.R <br>
Rscript paper_figures/figure3de.outlier.effect.size.R <br>
Rscript paper_figures/figure3.R

#### Figure 4
Rscript paper_figures/figure4a.uk10k.R <br>
Rscript paper_figures/figure4b.exac.enrichments.R <br>
Rscript paper_figures/figure4c.gene.list.enrichments.R <br>
Rscript paper_figures/figure4.R

#### Figure 5
Rscript getPosteriorsEval.R data/genomic_features.txt data/outliers.txt <br>
Rscript paper_figures/figure5b.R <br>
Rscript getPosteriorsApp.R data/genomic_features.txt data/outliers.txt data/postprobs_all.txt <br>
Rscript paper_figures/figure5c.R

## Supplemental figures

#### Number of rare variants per individual and PCA
You need to set the path to the downloaded expression covariates and subject annotations in the script below. <br> 
Rscript paper_figures/suppfig.number.rare.vars.pca.R

#### Improvement of replication of outliers across tissues by PEER correction
Rscript paper_figures/ExtendedDataFigure3.R

#### Distribution of the nubmer of genes with an multi-tissue outlier
You need to set the path to the subject annotations in the script below. <br> 
Rscript paper_figures/suppfig.number.outliers.per.individual.R

#### Single-tissue replication analysis controlling for sampling differences
Generated when running paper_figures/figure1b.outlier.sharing.R above

#### Overlap between single and multi-tissue outliers
Rscript paper_figures/suppfig.compare.single.multi.R

#### Single-tissue outlier rare variant enrichments
Generated when running paper_figures/figure2b.threshold.enrichments.R above

#### Multi-tissue outlier rare variant enrichments when excluding coding regions
Generated when running paper_figures/figure2a.count.enrichments.R above

#### Enrichment of functional genomic annotations among an expanded set of multi-tissue outliers
Rscript paper_figures/ExtendedDataFigure9.R

#### Association between ASE and RIVER scores
Rscript paper_figures/ExtendedDataFigure11.R

#### Correlation between test posterior probabilities and the fraction of tissues
Rscript paper_figures/main_RIVER_10CV.R <br>
Rscript paper_figures/ExtendedDataFigure12.R

#### Distribution of predictive scores for pathogenic variants and all variants
Rscript paper_figures/ExtendedDataFigure13.R

#### Expression levels for genes proximal to pathogenic variants
Rscript paper_figures/ExtendedDataFigure14.rpkm.R <br>
Rscript paper_figures/ExtendedDataFigure14.Zscores.R
