#!/usr/bin/bash

PEER_DIR=${1}
PEERs=${2}

# Function to compute residuals for given tissue
residCalculator(){
	line=${1}
	RPKM=${PEER_DIR}/${line}.rpkm.log2.ztrans.txt
	Covariates=${RAREVARDIR}/data/covariates/${line}_Analysis.covariates.txt
	Peer=${PEER_DIR}/${line}_Factors*/factors.tsv
	OUT=${PEER_DIR}/${line}.peer.top${PEERs}.ztrans.txt
	R -f PEER/calc_residuals_covs_peerless.R --slave --vanilla --args ${RPKM} ${Covariates} ${Peer} ${PEERs} ${OUT} &>${RAREVARDIR}/logs/resid.calculator.${line}.top${PEERs}.log
}

# Run above function for each tissue in parallel
# Requesting 10 cores
export PEER_DIR
export PEERs
export -f residCalculator
for line in `cat preprocessing/PEER/gtex_eqtl_tissues.txt`; do echo $line; done | parallel --jobs 10 residCalculator {1} 
