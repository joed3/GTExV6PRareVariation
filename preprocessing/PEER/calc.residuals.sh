#!/usr/bin/bash

set -o nounset -o errexit -o pipefail

PEER_DIR=${RAREVARDIR}/preprocessing/PEER

# Function to compute residuals for given tissue
residCalculator(){
	line=${1}
	RPKM=${PEER_DIR}/${line}.rpkm.log2.ztrans.txt
	Covariates=${RAREVARDIR}/data/eqtl_data/eQTLInputFiles/covariates/${line}_Analysis.covariates.txt
	Peer=${PEER_DIR}/${line}_Factors*/factors.tsv
	OUT=${PEER_DIR}/${line}.peer.ztrans.txt
	R -f ${PEER_DIR}/calc.residuals.core.R --slave --vanilla --args ${RPKM} ${Covariates} ${Peer} ${OUT} &>${RAREVARDIR}/logs/calc.residuals.${line}.log
}

# Run above function for each tissue in parallel
# Requesting 10 cores
export PEER_DIR
export -f residCalculator
for line in `cat ${PEER_DIR}/gtex_eqtl_tissues.txt`; do echo $line; done | parallel --jobs 10 residCalculator {1}
