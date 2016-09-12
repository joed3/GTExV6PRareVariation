#!/usr/bin/bash

PEER_DIR=${RAREVARDIR}/preprocessing/PEER

# Function of calculate PEER factors for a given tissue
factorCalculator(){
	line=${1}
	bash ${PEER_DIR}/calc.PEER.factors.single.tissue.sh &>${RAREVARDIR}/logs/calc.PEER.factors.${line}.log
}

# Calculate PEER factors for each tissue
# Run in parallel requesting 10 cores
export PEER_DIR
export -f factorCalculator
for line in $(cat ${PEER_DIR}/gtex_eqtl_tissues.txt); do echo $line; done | parallel --jobs 10 factorCalculator {1}
