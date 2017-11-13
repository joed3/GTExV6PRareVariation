# ===================================================================================================================
# This is a main script for running RIVER with varying thresholds in EDF 9a
#
# input:
#       1. genomic_features_thrdx.txt: a [N x (P+2)] matrix with values of genomic features (x: 1,1p5,2,2p5)
#          N is the total number of instances (individual and gene pairs). 
#          P is the total number of genomic features.
#          First two columns consist of individual IDs and gene names.
#       2. Zscores_thrdx.txt: a [N x (T+2)] matrix with Z-score values  (x: 1,1p5,2,2p5)
#          N is the total number of instances (individual and gene pairs). 
#          T is the total number of tissues considered
#          First two columns consist of individual IDs and gene names.
#       3. N2pairs.txt
#          a list of N2 pairs for evaluating RIVER and GAM (row: individual and gene instances, column1: individual ID, column2: gene name)
#
# output:
#       1. RIVER_thrdx.RData (x: 1,1p5,2,2p5)
#          store all data used for running RIVER with outliers based on x threshold 
#       1. RIVER_VaryingThrds.out.RData          
#          store all data used for this script
#
# ===================================================================================================================

# Load required packages
rm(list = ls())
library(pROC)

# Master directory
dir = Sys.getenv('RAREVARDIR')

# Recall required functions
source(paste(dir,'/paper_figures/RIVER.R',sep=""))
library(pROC)
library(glmnet)
library(dplyr)
library(data.table)

thrds_table = data.frame(values=c(1,1.5,2,2.5),labels=c("1","1p5","2","2p5"))
auc_table = data.frame(RIVER=matrix(NA,nrow(thrds_table),1), 
                       GAM=matrix(NA,nrow(thrds_table),1),
                       nlog10pval=matrix(NA,nrow(thrds_table),1))
sig_table = data.frame(RIVER=matrix(NA,nrow(thrds_table),1), 
                       GAM=matrix(NA,nrow(thrds_table),1))

corr_table = data.frame(Rcont=matrix(NA,nrow(thrds_table),1),pvalcont=matrix(NA,nrow(thrds_table),1),
                        Rdisc=matrix(NA,nrow(thrds_table),1),pvaldisc=matrix(NA,nrow(thrds_table),1),
                        npairs=matrix(NA,nrow(thrds_table),1))

for (level in 1:nrow(thrds_table)) {
  ## load genomic features and outlier status used for training RIVER
  G = as.data.frame(fread(paste(dir,"/data/genomic_features_thrd",thrds_table[level,"labels"],".txt",sep=""), 
                          sep='\t', header = TRUE, na.strings = "NaN")) # genomic features
  Gtemp = data.frame(G %>% mutate(rnum=row_number()))
  
  E = as.data.frame(fread(paste(dir,"/data/Zscores_thrd",thrds_table[level,"labels"],".txt",sep=""),
                          sep='\t', header = TRUE, na.strings = "NaN")) # z-scores
  
  ## load information of all rare variants (only 10K from TSS)
  GAll = as.data.frame(fread(paste(dir,"/data/gtex_AllRareSNVs.uniq.complete.txt",sep=""), 
                             sep='\t', header = TRUE, na.strings = "NaN")) 
  GAll = GAll[,c(1,2,3,4,5,7)] # indiv | gene | chr | pos | allele | nrv
  
  Gtemp = data.frame(Gtemp, exid=paste(Gtemp[,"gene"],";",Gtemp[,"indiv"],sep=""))
  
  GAll = data.frame(GAll,exid=paste(GAll[,"gene"],";",GAll[,"indiv"],sep=""))
  GAll[,c("indiv","gene")] <- NULL # remove indiv and gene since they are redundant 
  
  Gn = left_join(Gtemp, GAll, by="exid")
  
  ## import N2 pairs
  dups <- load(paste(dir,'/data/N2pairs.txt',sep=""))
  corr_table[level,"npairs"] = dim(dups)[1]
  
  # Compute |Median Z-scores|
  E_vals = apply(E[,3:(ncol(E)-1)], 1, function(x) abs(median(x, na.rm=TRUE)))
  
  # Define multi-tissue outliers
  # E_disc = as.numeric(E_vals >= 1.5)
  E_disc = as.numeric(E_vals >= as.numeric(thrds_table[level,"values"]))
  
  cat("\nThreshold = ", as.numeric(thrds_table[level,"values"]), "\n", sep="")
  cat("# of samples = ", dim(G)[1], "\n", sep="")
  cat("# of genomic features = ", dim(G)[2]-2, "\n", sep="")
  cat("% of outliers = ", round(sum(E_disc)/length(E_disc)[1],3), "\n", sep="")
  cat("# of N2 pairs = ", dim(dups)[1], "\n", sep="")
  
  corr_cont = cor.test(E_vals[dups[,1]],E_vals[dups[,2]],alternative="g",method="kendall")
  corr_table[level,"Rcont"] = corr_cont$estimate
  corr_table[level,"pvalcont"] = corr_cont$p.value
  
  cat("cont R: ", corr_cont$estimate, ", pval = ", corr_cont$p.value, "\n", sep ="")
  
  corr_disc = cor.test(E_disc[dups[,1]],E_disc[dups[,2]],alternative="g",method="kendall")
  corr_table[level,"Rdisc"] = corr_disc$estimate
  corr_table[level,"pvaldisc"] = corr_disc$p.value
  
  cat("disc R: ", corr_disc$estimate, ", pval = ", corr_disc$p.value, "\n", sep ="")
  
  basic_data = cbind(E_disc, G)
  rp = sample.int(dim(basic_data)[1])
  
  # Split data into trainig and test dataset
  train_inds = setdiff(rp, union(dups[,1], dups[,2]))
  test_inds = union(dups[,1], dups[,2])
  
  # Standardize genomic features of training data
  g_all = as.matrix(G[,3:ncol(G)])
  mean_col_g = apply(g_all,2,mean)
  sd_col_g = apply(g_all,2,sd)
  g_all <- (g_all-matrix(rep(mean_col_g,dim(g_all)[1]),byrow=TRUE,nrow=dim(g_all)[1]))/
    matrix(rep(sd_col_g,dim(g_all)[1]),byrow=TRUE,nrow=dim(g_all)[1])
  
  # Generate training data by holding out N2 pairs
  g_trng = as.matrix(G[train_inds,3:ncol(G)])
  g_trng <- (g_trng-matrix(rep(mean_col_g,dim(g_trng)[1]),byrow=TRUE,nrow=dim(g_trng)[1]))/
    matrix(rep(sd_col_g,dim(g_trng)[1]),byrow=TRUE,nrow=dim(g_trng)[1])
  
  # Search a best lambda from a multivariate logistic regression with outlier status with 10 cross-validation
  costs = c(100, 10, 1, .1, .01, 1e-3, 1e-4) # a list of candidate lambdas for L2 penalty in multivariate logistic regression
  cv.ll = cv.glmnet(g_trng, as.vector(E_disc[train_inds]), lambda=costs, 
                    family="binomial", alpha = 0, nfolds=10) # genomeic annotation model
  
  # Compute a P(FR | G) for all data
  pp = predict(cv.ll, g_all, s="lambda.min", type='response')
  
  # Train RIVER on training data
  theta = matrix(c(.99, .01, .3, .7), nrow=2, ncol=2) # initial theta
  pseudocount = 50 # pseudo count
  em.res <- integratedEM(g_trng, E_disc[train_inds], cv.ll$lambda.min, 
                         pseudocount, theta, cv.ll$glmnet.fit, costs)
  
  # Generate G data for test data (individual 1 from N2 pairs)
  g_test = as.matrix(G[dups[,2], 3:ncol(G)])
  
  # Generate test data (individual 1 from N2 pairs)
  g_test <- (g_test-matrix(rep(mean_col_g,dim(g_test)[1]),byrow=TRUE,nrow=dim(g_test)[1]))/
    matrix(rep(sd_col_g,dim(g_test)[1]),byrow=TRUE,nrow=dim(g_test)[1])
  
  # Compute P(FR | G, E)
  dup.post = testPosteriors(g_test, E_disc[dups[,2]], em.res)
  
  # Check performance of models with N2 pairs
  int.em.roc = roc(basic_data$E_disc[dups[,1]], dup.post$posterior[,2]) # RIVER
  sig_table[level,"RIVER"] = wilcox.test(dup.post$posterior[which(basic_data$E_disc[dups[,1]]==1),2],
                                         dup.post$posterior[which(basic_data$E_disc[dups[,1]]==0),2],
                                         alternative="greater")$p.value
  g.roc = roc(basic_data$E_disc[dups[,1]], pp[dups[,1]]) # genomic annotation model
  temp_pp = pp[dups[,1]]
  sig_table[level,"GAM"] = wilcox.test(temp_pp[which(basic_data$E_disc[dups[,1]]==1)],
                                       temp_pp[which(basic_data$E_disc[dups[,1]]==0)],
                                       alternative="greater")$p.value
  
  # roc.test(int.em.roc, g.roc)$p.value
  auc_table[level,"RIVER"] = int.em.roc$auc
  auc_table[level,"GAM"] = g.roc$auc
  auc_table[level,"nlog10pval"] = roc.test(int.em.roc, g.roc)$p.value
  
  cat('*** AUC (GAM - genomic annotation model): ',round(g.roc$auc,3),'\n')
  cat('    AUC (RIVER): ',round(int.em.roc$auc,3),'\n')
  cat('    P-value: ',format.pval(roc.test(int.em.roc, g.roc)$p.value,digits=4,
                                  eps=0.00001),'***\n')
  # =========== Save data
  save.image(file = paste(dir,"/data/RIVER_thrd",
                          as.character(thrds_table[level,"labels"]), ".out.RData",sep=""))
}
save.image(file = paste(dir,"/data/RIVER_VaryingThrds.out.RData",sep=""))
