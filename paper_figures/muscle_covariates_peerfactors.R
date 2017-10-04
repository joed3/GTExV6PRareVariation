# =======================================================================================================
# This is a script of 1st step for generating EDF1a.
#
# input:
#       1. subject_covariates.RData
#          subject covariates, each data matrix consists of N X P where N is the number of subjects and P is the number of covariates
#       2. sample_covariates.RData
#          sample covariates, each data matrix consists of N X P where N is the number of subjects and P is the number of covariates
#
# output:
#       1. muscle_subject_covariates_peerFactors.RData
#          adjusted R-squared values between PEER factors and each of subject covariates in muscle skeletal
#       2. muscle_samples_covariates_peerFactors.RData
#          adjusted R-squared values between PEER factors and each of sample covariates in muscle skeletal
#
# =======================================================================================================

# !/usr/bin/env Rscript

# Recall required packages
rm(list=ls(all=TRUE))

# Master directory
dir = Sys.getenv('RAREVARDIR')

# Recall required functions
library(parallel)

load("/data/subject_covariates.RData")
load("/data/sample_covariates.RData")
peer.fn <- dir("/data/peer/", full.names = T) # estimated peer factors for each of gene expression across 44 tissues separately

dat.peer <- lapply(peer.fn, function(x){
	d <- read.delim(x, header = T, row.names = 1)
	keep.rows <- grep("InferredCov", rownames(d))
	d <- t(d[keep.rows,])
	d
	})


dat.subject <- lapply(dat.peer, function(q,p){
	p <-p[gsub('[.]','-',rownames(q)),]
	## again 0 for low values
	p <- p[,which(sapply(p, function(y) sum(table(unique(y)))>1))]
	p <- p[,which(sapply(p, function(z) length(which(!is.na(z)))>=50))]
	p
	}, dat.subject)

dat.subject <- lapply(dat.subject, function(p){
	sample.names <- rownames(p)
	p <- as.data.frame(lapply(p, function(x){
		sapply(x, function(y){
			if(class(y) == "character"){
				y <- as.factor(y)
			}
			y
		})
		}))
	rownames(p) <- sample.names
	p
})


pve.function <- function(x,y){
	print("begin")
	y <- y[gsub('[.]','-',rownames(x)),]
	pve.vector <- matrix(NA, nrow = ncol(x), ncol = ncol(y))
	for(i in 1:ncol(x)){
		for(j in 1:ncol(y)){
			pve.vector[i,j] <- summary(lm(x[,i]~y[,j]))$adj.r.squared
		}
	}
	colnames(pve.vector) <- colnames(y)
	rownames(pve.vector) <- colnames(x)
	pve.vector
}
	
pve.explained <- pve.function(dat.peer[[29]], dat.subject[[29]])
#names(pve.explained) <- names(dat.samples)
save(pve.explained, file = paste("/data/muscle_subject_covariates_peerFactors.RData",sep=""))
rm(pve.explained)
pve.explained <- pve.function(dat.peer[[29]], dat.samples[[29]])
save(pve.explained, file = paste("/data/muscle_samples_covariates_peerFactors.RData",sep=""))
