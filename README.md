# GTExV6PRareVariation
Repository to reproduce analyses from the GTEx V6P Rare Variation Manuscript

# To run the code
## Download required files
Download from \<website\>: <br>
* processed data directory (referred to below as \<processed_data\>)

Download from http://www.gtexportal.org/home/datasets: <br>
* gencode.v19.genes.v6p_model.patched_contigs.gtf.gz
* GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz (gunzip it) <br>
* GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz (gunzip it) <br>
* GTEx_Data_V6_Annotations_SampleAttributesDS.txt

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

Other file: <br>
(Obtained directly from the Hall Lab.)
* gtex.lumpy.gs.svscore.low_conf.vcf.gz

## Setup
The code relies on an assumed directory structure.
Everything is under an upper level directory.
If you are using the processed files available online, set this upper-level directory to be the path to \<processed_data\>. <br>
Otherwise, set this to be any path where you want to scripts to output.

Either run or add the following to your .bashrc (omit trailing slashes for all directories): <br>
export RAREVARDIR=/<the path to the upper-level directory/>
export CADD_DIR=/<the path to the CADD files/>
export KG_DIR=/<the path to the 1000 genomes files/>
export ENCODE_MOTIF_DIR=/<the path to />/matches.txt.gz
export ER_DIR=/<the path to the Epigenomics Roadmap directory (sub directories are prom, enh, and dyadic)/>

If you are going to run everything from scratch, first make the directories that are assumed to exist. <br>
You can skip this step if you will use the processed data files downloaded from <\website\>.
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
Some of the steps reproduce processed files available on \<website\> and are marked as such.

# Pipeline
## Expression data correction and normalization
Generates processed data that can be downloaded from \<website\>. <br>

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
    --SAMPLE ${RAREVARDIR}/preprocessing/gtex_2015--1-12_samples_tissues.txt \ <br>
    --OUT ${RAREVARDIR}/preprocessing/PEER \ <br>
    --END .reads.txt

#### Make list of individuals and tissues
bash preprocessing/get_tissue_by_individual.sh

#### Run bash scripts to generate PEER corrected data (includes non-EAs) with covariates removed
bash PEER/goats_peer_pipeline_covs.sh ${RAREVARDIR}/preprocessing/PEER

#### Make flat files from raw RPKMs and PEER-corrected data for all tissues and individuals
(Takes 10 minutes or so to run) <br>
python preprocessing/gather_filter_normalized_expression.py <br>
python preprocessing/gather_filter_rpkm.py

## Preparing reference files used later.
Generates processed data that can be downloaded from \<website\>. <br>

Copy or move the GTEx annotation file (gencode.v19.genes.v6p_model.patched_contigs.gtf.gz) to ${RAREVARDIR}/reference.

bash preprocessing/process.reference.files.sh \<path to\>/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz \<path to\>/gencode.v19.annotation.gtf <br>
(relies on pad.gtf.exons.py, gtf2TSS.sh, and gtf2genebed.sh in preprocessing/)

## Outlier calling
Generates processed data that can be downloaded from \<website\>. <br>

#### Call multi-tissue outliers
(Uses multiple cores. Can set the number at the top of the script. Takes several minutes to run.) <br>
Rscript call_outliers/call_outliers_medz.R

#### Call single-tissue outliers
(Also takes several minutes to run.) <br>
python call_outliers/call_outliers_single_tissue.py

#### Compare single-tissue and multi-tissue outliers as well as get stats on each.
Rscript call_outliers/compare_single_multi_outliers.R

#### Run replication for single-tissue and multi-tissue outliers
(Both take a while to run. The multi-tissue replication uses multiple cores.)
Rscript multi_tissue_replication.R
Rscript single_tissue_replication.R

## Feature generation

#### Processing VCFs into bed files for each individual
bash vcf2bedfiles.sh \<path to \>/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz 
\<path to\>/gtex.lumpy.gs.svscore.low_conf.vcf.gz <br>
(relies on : <br>
* vcf2bedfiles_helper_processVCF.sh
* vcf2bedfiles_helper_processVCF_SV.sh
* vcf2bedfiles_helper_processVCFtoolsOutput.sh
* vcf2bedfiles_helper_processVCFtoolsOutput_CNV.sh
* compileCADDscores.sh
* extractCADDscores_ekt.py
*)

#### Extract features to be combined with the individual bed files
bash extract.1kg.AF.sh <br>
(relies on process.1kg.AF.py)

bash subset.CADD.features.sh <br>
(takes many hours)

bash TFBS_pipeline.sh

bash ER_pipeline.sh

## Main figures
#### Figure 1
bash pick.cartoon.example.sh <br>
Rscript figure1a.plot.cartoon.example.R <br>
Rscript figure1b.outlier.sharing.R <br>
Rscript figure1c.replication.rate.consistent.R

#### Figure 2

#### Figure 3

#### Figure 4

## Supplemental figures

#### Overlap between single and multi-tissue outliers
Rscript suppfig.compare.single.multi.R