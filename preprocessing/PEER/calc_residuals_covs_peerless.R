args <- commandArgs(trailingOnly=T)
if (length(args) != 5) {
  cat("Usage: R -f calc_residuals_covs_peerless.R --slave --vanilla --args RPKM COV PEER TOP.PEERs OUT\n", file=stderr())
  quit(status=2)
}


##Load require packages
require(data.table)
require(ggplot2)
require(reshape2)


#-------------- FUNCTIONS

#-------------- MAIN

# For testing
#expr_file = '/users/joed3/goats/data/gtex/rna-seq/PEER_genPCs_sex_wAAs/Adipose_Subcutaneous.rpkm.log2.ztrans.txt'
#peer_file = '/users/joed3/goats/data/gtex/rna-seq/PEER_genPCs_sex_wAAs/Adipose_Subcutaneous_Factors35/factors.tsv'
#covs_file = '/users/joed3/goats/data/gtex/eqtl_data/eQTLInputFiles/covariates/Adipose_Subcutaneous_Analysis.covariates.txt'

# Define arguments
expr_file = args[1]
covs_file = args[2]
peer_file = args[3]
top_peer_factors_to_remove = as.numeric(args[4])
out_file = args[5]


# Read in expression matrix
expr = read.table(expr_file, header = T, sep = '\t', row.names = 1)

# Read in covariates matrix
covs = read.table(covs_file, header = T, sep = '\t', row.names = 1)

# Replace periods in colnames with dashes
colnames(covs) = gsub('\\.', '-', colnames(covs))

# Reorder columns in covariates file to match expression matrix rows
covs = covs[, rownames(expr)]

# Remove platform covariate from analysis
covs = covs[-nrow(covs), ]

# Remove old PEER factors
old_indices = grep('InferredCov', rownames(covs))
covs = covs[-old_indices, ]

# Read in PEER factors
peer = read.table(peer_file, header = T, sep = '\t', row.names = 1)

# Replace periods in colnames with dashes
colnames(peer) = gsub('\\.', '-', colnames(peer))

# Reorder peer columns to match expression rows
# Only take the top N PEER factors
# If no PEER factors are to be removed, just remove the known covariates
if(top_peer_factors_to_remove > 0){
	peer = peer[1:top_peer_factors_to_remove, rownames(expr)]
	covs = rbind(covs, peer)
}

# Remove individuals with missing covariates
inds_to_keep = colSums(is.na(covs)) == 0
covs = covs[, inds_to_keep]
expr = expr[inds_to_keep, ]

# For each gene in the expression file, perform a linear regression 
# Keep residuals
resids = matrix(, ncol = ncol(expr), nrow = nrow(expr))
rownames(resids) = rownames(expr)
colnames(resids) = colnames(expr)

for(i in 1:ncol(expr)){
	#print(i)
	data = as.data.frame(cbind(expr[, i], t(covs)))
	colnames(data) = c('RPKM', rownames(covs))
	model = lm(RPKM ~ ., data = data)
	resids[, i] = model$residuals
}

# Center and scale, then transpose
resids = t(scale(resids))

# Write out the residuals
write.table(matrix(c('Id', colnames(resids)), nrow = 1), out_file, quote = F, row.names = F, col.names = F, sep = '\t')
write.table(resids, out_file, row.names = T, col.names = F, quote = F, sep = '\t', append = T)
