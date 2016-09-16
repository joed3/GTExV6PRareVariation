#!/usr/bin/Rscript

# R script to process RPKM and read count matrices so that columns match covariate files
# Prepare matrices for input to PEER

# Load require packages
require(data.table)
require(ggplot2)
require(reshape2)


#------------- FUNCTIONS

#------------- MAIN

dir = Sys.getenv('RAREVARDIR')
peer.dir = paste0(Sys.getenv('PEER_DIR'), '/')
covs.dir = paste0(dir, '/data/covariates/')

# Read in list of tissues with eQTL data
tissues = read.table('preprocessing/PEER/gtex_eqtl_tissues.txt', header = F, stringsAsFactors = F)[, 1]

# For each tissue, read in the RPKM and covariate files
# Subset and reorder columns of RPKM file to match covariate file
# Filter for genes with >= 10 individuals with RPKM > 0.1 and read count > 6
# Log2 + 2 transform the data, then z-transform
ind_filt = 10
rpkm_filt = 0.1
read_filt = 6

for(i in 1:length(tissues)){
	tissue = tissues[i]
	print(i)
	print(tissue)
	rpkm = as.data.frame(fread(paste0(peer.dir, tissue, '.rpkm.txt'), header = T))
	rownames(rpkm) = rpkm[, 1]
	rpkm = rpkm[, -1]
	reads = as.data.frame(fread(paste0(peer.dir, tissue, '.reads.txt'), header = T))
	rownames(reads) = reads[, 1]
	reads = reads[, -1]
	reads = reads[rownames(rpkm), names(rpkm)]
	covariates = read.table(paste0(covs.dir, tissue, '_Analysis.covariates.txt'), header = T, stringsAsFactors = F, row.names = 1)
	colnames(covariates) = gsub('\\.', '-', colnames(covariates))

	rpkm = rpkm[, colnames(covariates)]
	reads = reads[, colnames(covariates)]
	indices_to_keep = rowSums(rpkm > rpkm_filt & reads > read_filt) >= ind_filt
	rpkm = rpkm[indices_to_keep, ]
	rpkm_out = scale(t(log2(rpkm + 2)))
	write.table(rpkm_out, paste0(peer.dir, tissue, '.rpkm.log2.ztrans.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
}

