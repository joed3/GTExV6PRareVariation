# =======================================================================================================
# This is a script of 3rd step for generating EDF1b.
#
# input:
#       1. sample_covariates.RData
#          sample covariates, each data matrix consists of N X P where N is the number of subjects and P is the number of covariates
#       2. subject_covariates.RData
#          subject covariates, each data matrix consists of N X P where N is the number of subjects and P is the number of covariates
#
# output:
#       1. superheat.sample.RData and superheat.subject.RData
#          adjusted R-qaured values of either sample and subject covariates in each regressed-out data with PEER
#
# =======================================================================================================

# !/usr/bin/env Rscript

# Recall required packages
rm(list=ls(all=TRUE))

# Master directory
dir = Sys.getenv('RAREVARDIR')

# Recall required functions
library(reshape2)
library(gtools)

load("/data/sample_covariates.RData")
fn <- mixedsort(dir(path = "/data/", pattern = "pve_sample", full = T))
pve.mat <- lapply(fn, function(x) {
  load(x)
  pve.explained
})

names(pve.mat) <- names(dat.samples)
combined.pve <- plyr::ldply(pve.mat,rbind)
rownames(combined.pve) <- combined.pve$.id
combined.pve <- abs(combined.pve[,-1])
d <- scale(t(combined.pve))
ord <- hclust( dist(d, method = "euclidean"), method = "average" )$order
ord

col.heat <- rev(c("#561a65ff", "#8997b9ff","#90c7c5ff",
                  "#90c7c5ff", "#85d689ff", "#bee9c0ff", "#fdee72ff"))

tmp.super <- t(combined.pve[,ord])
save(tmp.super, file = "/data/superheat.sample.RData")

################## subject specifc

load("/data/sample_covariates.RData")
fn <- mixedsort(dir(path = "/data/", pattern = "pve_subject", full = T))

pve.mat <- lapply(fn, function(x) {load(x)
  pve.explained
})

combined.pve <- plyr::ldply(pve.mat,rbind)


rownames(combined.pve) <- names(dat.samples)
#combined.pve <- abs(combined.pve[,-1])

d <- scale(t(combined.pve))
ord <- hclust( dist(d, method = "euclidean"), method = "average" )$order
ord

col.heat <- rev(c("#561a65ff", "#8997b9ff","#90c7c5ff",
                  "#90c7c5ff", "#85d689ff", "#bee9c0ff", "#fdee72ff"))

tmp.super <- t(combined.pve[,ord])
save(tmp.super, file = "/data/superheat.subject.RData")
