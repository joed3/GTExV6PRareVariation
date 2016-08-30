#!/usr/bin/env Rscript

### Master directory
dir = Sys.getenv('RAREVARDIR')

# Load required packages
require(ggplot2)
require(reshape2)
require(data.table)
require(plyr)
require(doMC)

###
### Setup parallel processing
###
doMC::registerDoMC(cores=12)

#------------- FUNCTIONS

###
### Functions for calculating number of tissues/sample and meta-analysis z-score
###
meta.n = function(values) {
	length(values) - sum(is.na(values))
}

meta.median = function(values) {
	median(values, na.rm=T)
}

meta.analysis = function(x) {
	samples=colnames(x)[3:ncol(x)]
	y = t(x[,3:ncol(x)]) # individuals (rows) x tissues (columns)
	n = apply(y, 1, meta.n)
	m1 = apply(y, 1, meta.median)
	data.frame(sample = samples, n.tissues = n, median.z = m1)
}

# Function tp pick most extreme outlier per gene
meta.extreme = function(x, tissue.filter = 5) {
	max.index = which(abs(x$median.z) == max(abs(x$median.z[x$n.tissues >= tissue.filter]), na.rm = T))
	x[max.index, ]
}

# Function to pick random outlier per gene
meta.random = function(x, tissue.filter = 5){
	index = sample(which(x$n.tissues >= tissue.filter), 1)
	x[index, ]
}

#------------- MAIN

# Load flat file with filtered and normalized expression data
data = fread(paste0(dir, '/preprocessing/gtex_2015-01-12_normalized_expression.txt'), header = T)

setkey(data, Gene)

# Get list of individuals with data
individs = colnames(data)[3:ncol(data)]

# Get list of tissues
tissues = unique(sort(data$Tissue))

# Read in list of GENCODE genes with types
genes.types = read.table(paste(dir, '/reference/gencode.v19.genes.v6p.patched_contigs_genetypes_autosomal.txt'), sep = '\t', header = F, stringsAsFactors = F)

# Filter for protein_coding and lincRNA genes
types.to.keep = c('protein_coding', 'lincRNA')
genes.types = genes.types[genes.types[, 2] %in% types.to.keep, ]
data = data[data$Gene %in% genes.types[, 1], ]

# Define tissue sizes to test replication
tissue.sizes = seq(10, 30, by = 5)

# Define tissues to hold out for replication (a consistent set for all discovery sizes)
set.seed(10)
replication.size = 10
replication.tissues = sample(tissues, replication.size)
replication.data.full = data[data$Tissue %in% replication.tissues, ]
discovery.tissues.full = tissues[!(tissues %in% replication.tissues)]

# For each tissue size, sample discovery set of that size
# Use remaining tissues for replication set
# Repeat at each sample size N times
discovery.thresh = 2
replication.thresholds = c(1, 2)
n.samps = 10
numerators.real = matrix(, ncol = n.samps * length(replication.thresholds), nrow = length(tissue.sizes))
denominators.real = matrix(, ncol = n.samps * length(replication.thresholds), nrow = length(tissue.sizes))
numerators.back = matrix(, ncol = n.samps * length(replication.thresholds), nrow = length(tissue.sizes))
denominators.back = matrix(, ncol = n.samps * length(replication.thresholds), nrow = length(tissue.sizes))
for(i in 1:length(tissue.sizes)){
	for(j in 1:n.samps){
		discovery.tissues = sample(discovery.tissues.full, tissue.sizes[i])
		discovery.data = data[data$Tissue %in% discovery.tissues, ]

		# Will the real outlier please stand up?
		# Discover outliers in discovery set
		discovery.results = ddply(discovery.data, .(Gene), meta.analysis, .parallel = TRUE)
		discovery.outliers = ddply(discovery.results, .(Gene), meta.extreme, .parallel = TRUE)
		discovery.outliers = discovery.outliers[abs(discovery.outliers$median.z) >= discovery.thresh, ]
		discovery.outliers$Pair = paste(discovery.outliers$Gene, discovery.outliers$sample, sep = '_')
		discovery.outliers$Outlier = (abs(discovery.outliers$median.z) >= discovery.thresh) + 0

		# Generate replication data
		replication.data = replication.data.full[replication.data.full$Gene %in% discovery.outliers$Gene, ]
		replication.results = ddply(replication.data, .(Gene), meta.analysis, .parallel = TRUE)
		replication.results = replication.results[replication.results$n.tissues > 0, ]
		replication.results$Pair = paste(replication.results$Gene, replication.results$sample, sep = '_')

		# Generate data for background replication rate
		discovery.results.back = discovery.results[discovery.results$Gene %in% discovery.outliers$Gene, ]
		discovery.outliers.back = ddply(discovery.results.back, .(Gene), meta.random, .parallel = TRUE)
		discovery.outliers.back$Pair = paste(discovery.outliers.back$Gene, discovery.outliers.back$sample, sep = '_')

		# Calculate replication rate at various thresholds 
		for(k in 1:length(replication.thresholds)){
			out.thresh = replication.thresholds[k]
			print(paste('Tissue size ', tissue.sizes[i], ', sample ', j, ', threshold ', out.thresh, sep = ''))
			# Calculate replication
			# For now, will use percentage of outliers with Median-Z greater than 2 in replication set
			replication.results$Outlier = (abs(replication.results$median.z) >= out.thresh) + 0
			replication.outliers = replication.results[replication.results$Pair %in% discovery.outliers$Pair, ]
			outliers.merged = merge(discovery.outliers, replication.outliers, by = 'Pair')
			numerators.real[i, j + (k - 1) * n.samps] = sum(outliers.merged$Outlier.y == 1 & sign(outliers.merged$median.z.x) == sign(outliers.merged$median.z.y)) 
			denominators.real[i, j + (k - 1) * n.samps] = nrow(outliers.merged)

			# What is the background replication rate?
			replication.outliers.back = replication.results[replication.results$Pair %in% discovery.outliers.back$Pair, ]
			outliers.merged.back = merge(discovery.outliers.back, replication.outliers.back, by = 'Pair')
			numerators.back[i, j + (k - 1) * n.samps] = sum(outliers.merged.back$Outlier == 1 & sign(outliers.merged.back$median.z.x) == sign(outliers.merged.back$median.z.y)) 
			denominators.back[i, j + (k - 1) * n.samps] = nrow(outliers.merged.back)
		}
	}
}

# Make the replication plot
rep.types = c('Observed', 'Background')

replication.real = data.frame(SIZE = rep(tissue.sizes, n.samps), 
	NUM = as.vector(numerators.real), 
	DEN = as.vector(denominators.real), 
	THRESH = factor(rep(replication.thresholds, each = n.samps * length(tissue.sizes)), levels = replication.thresholds))
replication.real$Type = factor(rep(rep.types[1], nrow(replication.real)), levels = rep.types)

replication.back = data.frame(SIZE = rep(tissue.sizes, n.samps), 
	NUM = as.vector(numerators.back), 
	DEN = as.vector(denominators.back), 
	THRESH = factor(rep(replication.thresholds, each = n.samps * length(tissue.sizes)), levels = replication.thresholds))
replication.back$Type = factor(rep(rep.types[2], nrow(replication.back)), levels = rep.types)

replication = rbind(replication.real, replication.back)

# Save workspace image 
save(replication, replication.thresholds, file = paste0(dir, '/data/figure1c.replication.rate.consistent.RData'))
