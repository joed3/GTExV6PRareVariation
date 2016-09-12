args <- commandArgs(trailingOnly=T)
if (length(args) < 4) {
  cat("Usage: R -f calc.PEER.factors.core.R --slave --vanilla --args TraitsFileName MaxFactorsN MaxIterations BoundTol VarTol e_pa e_pb a_pa a_pb OutDir\n", file=stderr())
  quit(status=2)
}

library(peer)

GetPeerModel <- function(m1, MaxFactorsN, MaxIterations, a_pa, a_pb, e_pa, e_pb){
	model = PEER()
	# set data and parameters
	cat("PEER_setNk\n")
	res <- PEER_setNk(model, MaxFactorsN)
	if (!is.null(res)) RiseError("not correct")
	cat("PEER_setPhenoMean\n")
	PEER_setPhenoMean(model, as.matrix(m1))
	if (!is.null(res)) RiseError("not correct")
	# set priors (these are the default settings of PEER)
	cat("PEER_setPriorAlpha\n")
	PEER_setPriorAlpha(model,a_pa,a_pb);
	if (!is.null(res)) RiseError("not correct")
	cat("PEER_setPriorEps\n")
	PEER_setPriorEps(model,e_pa, e_pb);
	if (!is.null(res)) RiseError("not correct")
	cat("PEER_setNmax_iterations\n")
	PEER_setNmax_iterations(model,MaxIterations)
	if (!is.null(res)) RiseError("not correct")
	return(model)
}

RiseError <- function(Statement){
	cat(sprintf("%s\n", Statement))
	quit(status=2)
}

# Number of factors does not have an effect on PEER, as the unused
# factors are switched off. 
FactorsNComparison <- function(m1, MaxIterations, a_pa, a_pb, e_pa, e_pb, Factors){
    pdf("FactorsVarianceAll.pdf",width=8,height=8)
    # plot factor weight variances for a large set of maximum number K of factors, and see if K has an effect
    count=0
    for (k in Factors){
    	cat(sprintf("Factors%s\n", k))
	count=count+1
        model = GetPeerModel(m1, k, MaxIterations, a_pa, a_pb, e_pa, e_pb)
        PEER_update(model)
        if (count==1) PlotType <- "p"
        else PlotType <- "l"
	plot(1.0 / PEER_getAlpha(model),xlab="Factors", ylab="Factor relevance", type=PlotType,xlim=c(0,length(PEER_getAlpha(model))), ylim=c(0,1000), main=sprintf("Factors: %s\n", k))
    }
    dev.off()
}

# Prior on the noise level or factor weight precisions affects factor inference
# In general the prior parameters for Alpha and Eps can be tweaked to explain some
# amount of variability in the data. In many settings, PEER is robust to changes in priors
# for several orders of magnitude
EpsPriorComparison <- function(m1, MaxIterations, a_pa, a_pb, MaxFactorsN){
    # plot factor weight variances for a large set of maximum number K of factors, and see if K has an effect
    for(e_pa in c(0.0001, 0.1, 1000)){
        for(e_pb in c(0.1,10,1000)){
            model = GetPeerModel(m1, MaxFactorsN, MaxIterations, a_pa, a_pb, e_pa, e_pb)
            PEER_setPriorEps(model,e_pa, e_pb);
            PEER_update(model)
            print(paste("Eps pa=", e_pa, "pb=", e_pb, "mean(residuals^2)=",mean((PEER_getResiduals(model))**2)))
        }
    }
}

WritePeerDataCore <- function(m, samples_id_v, transcript_id_v, out_file_name_tag){
	m <- data.frame(samples_id_v,m, stringsAsFactors=F)
	colnames(m) <- c("ID",transcript_id_v)
	#write.table(m, file=sprintf("%s.Transposed.tsv", out_file_name_tag), sep="\t", row.names=F, quote=F)
	m <- t(m)
	write.table(m, file=sprintf("%s.tsv", out_file_name_tag), sep="\t", col.names=F, quote=F)
}


WritePeerData <- function(m1, model, MaxIterations, a_pa, a_pb, e_pa, e_pb, MaxFactorsN){
	#factors
	X <- PEER_getX(model)
	#weights
	W <- PEER_getW(model)
	#ARD parameters, precision
	Alpha <- PEER_getAlpha(model)
	#get corrected dataset
	Yc <- PEER_getResiduals(model)

	SamplesIds <- row.names(m1)
	TranscriptIds <- colnames(m1)
	FactorsTags <- sprintf("Factors%s", 1:ncol(X))

	WritePeerDataCore(W, TranscriptIds, FactorsTags, "weigths")
	WritePeerDataCore(X, SamplesIds, FactorsTags, "factors")
	WritePeerDataCore(Alpha,FactorsTags,"precision", "precision")
	WritePeerDataCore(Yc, SamplesIds, TranscriptIds, "residuals")
	
	cat(sprintf("%s\t%s\n", "ID", "data"), file="PeerInfos.tsv")
	cat(sprintf("%s\t%s\n", "MaxIterations", MaxIterations), file="PeerInfos.tsv", append=T)
	cat(sprintf("%s\t%s\n", "MaxFactorsN", MaxFactorsN), file="PeerInfos.tsv", append=T)
	cat(sprintf("%s\t%s\n", "e_pa", e_pa), file="PeerInfos.tsv", append=T)
	cat(sprintf("%s\t%s\n", "e_pb", e_pb), file="PeerInfos.tsv", append=T)
	cat(sprintf("%s\t%s\n", "a_pa", a_pa), file="PeerInfos.tsv", append=T)
	cat(sprintf("%s\t%s\n", "a_pb", a_pb), file="PeerInfos.tsv", append=T)
	cat(sprintf("ResidualsSquareMean\t%s\n", mean((Yc)**2)), file="PeerInfos.tsv", append=T)
}

PlotFactorsVariance <- function(Alpha, tissue){
	# plot variance of factors - in this case, we expect
	# a natural elbow where there are 5 active factors, as 5 were simulated
	pdf(paste(tissue, "_FactorsVariance.pdf", sep = ''),width=8,height=8)
	plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance", main="")
	lines(1.0 / Alpha,xlab="Factors", ylab="Factor relevance", main="")
	dev.off()
}

PlotFactorsWeights <- function(model, tissue){
	pdf(paste(tissue, '_FactorsWeights.Scatterplot.pdf', sep = ''))
	PEER_plotModel(model)
	dev.off()
}



###################################
### setting parameters
###################################
#TraitsFileName <- "/SHARE/USERFS/els7/users/IRGB_CNR/Tmp/examples/data/expression.csv"
#MaxFactorsN <- 20
#MaxIterations <- 100
TraitsFileName <- args[1]
MaxFactorsN <- as.numeric(args[2])
MaxIterations <- as.numeric(args[3])


# PEER finishes if the increase in lower bound on the model evidence ceases
# to change, or the variance of the residuals has stabilised.
# Set  limiting values (tolerances)
# You can keep the bound tolerance fairly high, but should keep the variation
# tolerance quite low compared to the variance of the expression matrix.
# default
#BoundTol <- 0.001
#VarTol <- 0.00001
BoundTol <- as.numeric(args[4])
VarTol <- as.numeric(args[5])

# Prior parameters on the noise and weight precision distributions can also
# be changed. As these are both gamma distributed, you can specify
# the a and b parameters of both
### eps, noise
#e_pa <- 0.1
#e_pb <- 10.
#### alpha, weights
#a_pa <- 0.001
#a_pb <- 0.1

e_pa <- as.numeric(args[6])
e_pb <- as.numeric(args[7])
a_pa <- as.numeric(args[8])
a_pb <- as.numeric(args[9])
OutDir <- args[10]
tissue <- args[11]


cat(sprintf("TraitsFileName: %s\n", TraitsFileName))
cat(sprintf("MaxFactorsN: %s\n", MaxFactorsN))
cat(sprintf("MaxIterations: %s\n", MaxIterations))
cat(sprintf("BoundTol: %s\n", BoundTol))
cat(sprintf("VarTol: %s\n", VarTol))
cat(sprintf("e_pa: %s\n", e_pa))
cat(sprintf("e_pb: %s\n", e_pb))
cat(sprintf("a_pa: %s\n", a_pa))
cat(sprintf("a_pb: %s\n", a_pb))


###################################
### loading data
###################################
cat ("loading data\n")
m1 <- read.table(TraitsFileName, header=T, row.names=1, sep = '\t')
#m1 <- read.delim(TraitsFileName, header=F, sep=",")
cat(sprintf("loaded %s samples %s traits\n", dim(m1)[1], dim(m1)[2]))

dir.create(OutDir)
setwd(OutDir)

###################################
### Peer
##################################
model <- GetPeerModel(m1, MaxFactorsN, MaxIterations, a_pa,a_pb, e_pa, e_pb)
cat("PEER_update\n")
PEER_update(model)

###################################
### get results
###################################
cat("writing data\n")
WritePeerData(m1, model, MaxIterations, a_pa, a_pb, e_pa, e_pb, MaxFactorsN)
PlotFactorsVariance(PEER_getAlpha(model), tissue)
PlotFactorsWeights(model, tissue)
#FactorsNComparison(m1, MaxIterations, a_pa, a_pb, e_pa, e_pb, 1:10)
#EpsPriorComparison(m1, MaxIterations, a_pa, a_pb, MaxFactorsN)

cat("calc.PEER.factors.core.R DONE\n")


