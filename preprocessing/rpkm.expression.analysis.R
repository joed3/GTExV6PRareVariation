#!/usr/bin/env Rscript
# Load packages

require(reshape2)
require(plyr)
require(data.table)
require(doMC)
require(matrixStats)

# Register backend for parallelization
doMC::registerDoMC(cores = 10)

#------------ FUNCTIONS


#-------------- MAIN

# Read in list of individuals and tissues
samples = scan('preprocessing/gtex_2015-01-12_individuals_all_normalized_samples.txt', what = character())
tissues = scan('preprocessing/gtex_2015-01-12_tissues_all_normalized_samples.txt', what = character())

# Read in the flat file 
rpkm = fread('preprocessing/gtex_2015-01-12_rpkm.txt')

# Gene and Tissue IDs as keys
keys = c('Tissue', 'Gene')
setkeyv(rpkm, keys)

# Read in list of all expressed genes 
expressed.genes = read.table('preprocessing/gtex.expressed.genes.txt', sep = '\t', header = F, stringsAsFactors = F)[, 1]

# Filter RPKM flat file for only these expressed genes
rpkm = rpkm[Gene %in% expressed.genes]

# Calculate mean gene expression levels
# Takes a few minutes
rpkm_mean = ddply(rpkm, .(Gene, Tissue), function(x) mean(unlist(x[1, 3:ncol(x)]), na.rm = T), .parallel = T)

# Unmelt the datasets to yield matrices of genes by tissues
rpkm_mean_mat = dcast(rpkm_mean, Gene ~ Tissue)

rownames(rpkm_mean_mat) = rpkm_mean_mat$Gene

rpkm_mean_mat = rpkm_mean_mat[, -1]

rpkm_mean_mat = as.matrix(rpkm_mean_mat)

# Calculate mean/median/sd across tissues
tissue_means = rowMeans(rpkm_mean_mat)
tissue_sds = rowSds(rpkm_mean_mat)
tissue_meds = apply(rpkm_mean_mat, 1, median, na.rm = T)

tissue_stats = data.frame(GENE = rownames(rpkm_mean_mat), MEAN = tissue_means, MEDIAN = tissue_meds, SD = tissue_sds)

# Write out tissue stats and raw matrices
write.table(tissue_stats, 'data/genes.rpkm.summary.stats.txt', quote = F, sep = '\t', col.names = T, row.names = F)
write.table(rpkm_mean_mat, 'data/genes.rpkm.mean.mat.txt', quote = F, sep = '\t', col.names = T, row.names = T)

