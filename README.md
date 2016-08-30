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
gencode.v19.annotation.gtf.gz (gunzip it)

## Setup
The code relies on an assumed directory structure.
Everything is under an upper level directory.
If you are using the processed files available online, set this upper-level directory to be the path to \<processed_data\>. <br>
Otherwise, set this to be any path where you want to scripts to output.

Either run or add the following to your .bashrc: <br>
export RAREVARDIR=<the path to the upper-level directory without trailing slash>


If you are going to run everything from scratch, first make the directories that are assumed to exist. <br>
You can skip this step if you will use the processed data files downloaded from <\website\>.
mkdir ${RAREVARDIR}/preprocessing <br>
mkdir ${RAREVARDIR}/preprocessing/PEER <br>
mkdir ${RAREVARDIR}/reference <br>
mkdir ${RAREVARDIR}/data <br>
mkdir ${RAREVARDIR}/data/medz <br>
mkdir ${RAREVARDIR}/data/singlez <br>

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

bash preprocessing/process.reference.files.sh \<path to downloaded file\>/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz \<path to downloaded file\>/gencode.v19.annotation.gtf <br>
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





## Supplemental figures

#### Overlap between single and multi-tissue outliers
Rscript suppfig.compare.single.multi.R