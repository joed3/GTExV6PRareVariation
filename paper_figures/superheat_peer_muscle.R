# =======================================================================================================
# This is a script of 2nd step for generating EDF1a.
#
# input:
#       1. cov_for_heatmap.txt
#          a list of names for both sample and subject covariates
#
# output:
#       EDF1a figure
#
# Note that the final figure was slightly modified by using inkscape for visualization purposes
#
# =======================================================================================================

#!/usr/bin/env Rscript

# Recall required packages
rm(list=ls(all=TRUE))

# Master directory
dir = Sys.getenv('RAREVARDIR')

# Recall required functions
library(gplots)

inputargs <- commandArgs(TRUE)
inputfn <- inputargs[1]
savefn <- inputargs[2]
load(inputfn)

cov.names <- read.delim("/data/cov_for_heatmap.txt",sep = "\t", row.names =1)
tmp.super <- data.frame(pve.explained)
tmp.super[tmp.super < 0.0] <- 0.0
print(range(tmp.super, na.rm = T))
tmp.super <- t(tmp.super[,rev(names(sort(sapply(tmp.super, function(p) mean(p, na.rm = T)))))[1:20]])
print(rownames(tmp.super)[which(!rownames(tmp.super) %in% cov.names$V1)])
cov.names <- cov.names[which(cov.names$V1 %in% rownames(tmp.super)),]
rownames(cov.names) <- cov.names$V1
cov.names <- cov.names[rownames(tmp.super),]
rownames(tmp.super) <- cov.names$V2
colnames(tmp.super) <- c(1:ncol(tmp.super))

pdf(savefn, height = 13, width = 15, useDingbats = F)
heatmap.2(as.matrix(tmp.super), col = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(40),
	scale="none",
          margins=c(8,8), # ("margin.Y", "margin.X")
          trace='none', 
          symkey=FALSE, 
          symbreaks=FALSE, 
          dendrogram='none',
          density.info='none', 
          denscol="black",
	  breaks = seq(0, 0.7, length.out = 41),
         # cexRow = 1.4,
	  keysize=0.3, 
          key.xlab="",
	  #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(3,0,5,6.5)),
	  key.title = "",
	  # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(7,8,9),c(5,1,2), c(6, 4, 3)), lhei=c(0.5,3.5, 0.5), lwid=c(0.5, 6, 1)
)
dev.off()
