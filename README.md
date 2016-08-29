# GTExV6PRareVariation
Repository to reproduce analyses from the GTEx V6P Rare Variation Manuscript

# To run the code
## Download required files
Download from \<website\>: <br>
* processed_data directory

Download from http://www.gtexportal.org/home/datasets (and gunzip): <br>
* GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz <br>
* GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz <br>

First make the directories that are assumed to exist <br>
mkdir ../preprocessing <br>
mkdir ../preprocessing/PEER <br>

Then you can run the code below. <br>
Some of the steps reproduce processed files available on \<website\> and are marked as such.

# Pipeline
## Expression data correction and normalization
Generates processed data that can be downloaded from \<website\>
If you skip this step, copy/move the following files from \<path to processed_data\> to ../preprocessing:
* gtex_2015-01-12_rpkm.txt 
* gtex_2015-01-12_normalized_expression.txt

#### Generate rpkm and read count matrices from GTEx V6P combined file
python preprocessing/split_by_tissues.py \ <br>
    --GTEX \<path to count data\>/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct \ <br>
    --SAMPLE \<path to processed_data directory\>/gtex_samples_tissues.txt \ <br>
    --OUT ../preprocessing/PEER \ <br>
    --END .rpkm.txt <br><br>

python preprocessing/split_by_tissues.py \ <br>
    --GTEX \<path to count data\>/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct \ <br>
    --SAMPLE \<path to processed_data directory\>/gtex_samples_tissues.txt \ <br>
    --OUT ../preprocessing/PEER \ <br>
    --END .reads.txt

#### Make list of individuals and tissues
bash preprocessing/get_tissue_by_individual.sh

#### Run bash scripts to generate PEER corrected data (includes non-EAs) with covariates removed
bash PEER/goats_peer_pipeline_covs.sh ../preprocessing/PEER

#### Make flat files from raw RPKMs and PEER-corrected data for all tissues and individuals
(Takes 10 minutes or so to run) <br>
python preprocessing/gather_filter_normalized_expression.py <br>
python preprocessing/gather_filter_rpkm.py

