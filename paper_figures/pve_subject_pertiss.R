# =======================================================================================================
# This is a script of 2nd step for generating EDF1b.
#
# input:
#       1. expression_component_removed.Rdata
#          PEER-corrected expression data per tissue, each data matrix consists of N x G where N is the number of subjects and G is the number of genes available
#       2. sample_covariates.RData
#          sample covariates, each data matrix consists of N X P where N is the number of subjects and P is the number of covariates
#       3. subject_covariates.RData
#          subject covariates, each data matrix consists of N X P where N is the number of subjects and P is the number of covariates
#
# output:
#       1. pve_subjects_i.RData 
#          proportion of variation explained by each covariate in i tissue
#
# =======================================================================================================

# !/usr/bin/env Rscript

# Recall required packages
rm(list=ls(all=TRUE))

# Master directory
dir = Sys.getenv('RAREVARDIR')

# Recall required functions
library(parallel)

load("/data/expression_component_removed.Rdata")

#tiss <- as.numeric(commandArgs(TRUE))
#tiss <- as.numeric("8")

for (tiss in (1:44)) {
  load("/data/subject_covariates.RData")
  load("/data/sample_covariates.RData")
  
  dat.subject <- dat.subject[gsub('[.]','-',rownames(expr.removed[[tiss]])),]
  
  ## again 0 for low values
  dat.subject <- dat.subject[,which(sapply(dat.subject, function(y) sum(table(unique(y)))>1))]
  dat.subject <- dat.subject[,which(sapply(dat.subject, function(z) length(which(!is.na(z)))>=50))]
  
  dat.subject <- as.data.frame(lapply(dat.subject, function(x){
    sapply(x, function(y){
      if(class(y) == "character"){
        y <- as.factor(y)
      }
      y
    })
  }))
  
  rownames(dat.subject) <- rownames(expr.removed[[tiss]])
  
  pve.function <- function(x,y){
    print("begin")
    pve.vector <- vector(length = ncol(y))
    for(j in 1:ncol(y)){
      print(paste(j))
      select.obs <- which(!is.na(y[,j]))
      x.sub <- scale(x[select.obs,])
      y.sub <- y[select.obs,]
      tts = (norm(x.sub, type = "F"))^2
      res <- lm(x.sub~y.sub[,j])$residuals
      # computing residual sum of squares
      dof <- summary(lm(x.sub[,1]~y.sub[,j]))$df
      rss <- (norm(res, type = "F"))^2
      # sum of squares explained
      r.sq <- (tts - rss)/tts
      adj.r.sq <- r.sq - ((1-r.sq) * (dof[1]/(dof[2]-1)))
      pve.vector[j] <- adj.r.sq
      #pve.vector[j] <- r.sq
      
    }
    names(pve.vector) <- colnames(y)
    pve.vector
  }
  
  pve.explained <- pve.function(expr.removed[[tiss]], dat.subject)
  save(pve.explained, file = paste("/data/pve_subjects_", tiss,".Rdata", sep = ""))
}

