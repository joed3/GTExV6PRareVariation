#!/usr/bin/Rscript

# Load required packages 
require(reshape2)
require(plyr)
require(ggplot2)

#--------------- FUNCTIONS

#--------------- MAIN

dir = Sys.getenv('RAREVARDIR')

# Define minimum P-value
minP = 1e-313

# Read in the Metasoft BF corrected data
meta = read.table(paste0(dir, '/data/metasoft/gtex.metasoft.v6p.selected.txt'), sep = '\t', header = T, stringsAsFactors = F)

# Read in Metasoft tissue order
tissues = read.table(paste0(dir, '/data/metasoft/Metasoft_tissue_order.txt'), sep = '\t', header = F, stringsAsFactors = F)[, 1]

# Convert M-values to eQTL states
metasoft.thresh = 0.9
mvalues = meta[, 62:105]
evalues = (mvalues >= metasoft.thresh) + 0
ntissues.egenes = rowSums(evalues, na.rm = T)
ntissues.tested = rowSums(!is.na(evalues))

# Make data frame for plotting
egenes.df = data.frame(Gene = meta$GENE, Tested = ntissues.tested, eGenes = ntissues.egenes, Isq = meta$I_SQUARE, Q = meta$Q, 
	Pq = meta$PVALUE_Q, TauSq = meta$TAU_SQUARE, Pvalue = meta$PVALUE_RE2)

# Add RPKM data across tissues to the data frame
rpkm = read.table(paste0(dir, '/data/genes.rpkm.summary.stats.txt'), sep = '\t', header = T, stringsAsFactors = F)
egenes.df$RPKM = NA
for(i in 1:nrow(egenes.df)){
	if(egenes.df$Gene[i] %in% rpkm$GENE){
		egenes.df$RPKM[i] = rpkm$MEAN[rpkm$GENE == egenes.df$Gene[i]]
	}
}

# Write out egenes.df sorted by PVALUE_RE2
egenes.df = egenes.df[order(egenes.df$Pvalue), ]

write.table(egenes.df, paste0(dir, '/data/GTExReleaseV6PMetasoft.summary.txt'), sep = '\t', col.names = T, row.names = F, quote = F)




