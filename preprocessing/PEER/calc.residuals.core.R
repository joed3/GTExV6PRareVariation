args <- commandArgs(trailingOnly=T)
if (length(args) != 4) {
  cat("Usage: R -f calc.residuals.core.R --slave --vanilla --args RPKM COV PEER OUT\n", file=stderr())
  quit(status=2)
}


##Load require packages
require(data.table)
require(ggplot2)
require(reshape2)


#-------------- FUNCTIONS

#-------------- MAIN

# Define arguments
expr_file = args[1]
covs_file = args[2]
peer_file = args[3]
out_file = args[4]

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
peer = peer[, rownames(expr)]

# Combine covariates and PEER factors
covs = rbind(covs, peer)

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
