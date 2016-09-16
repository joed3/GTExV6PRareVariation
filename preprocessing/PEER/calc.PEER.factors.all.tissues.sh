#!/usr/bin/bash

set -o nounset -o errexit -o pipefail

# Function of calculate PEER factors for a given tissue
factorCalculator(){
	line=${1}
	bash preprocessing/PEER/calc.PEER.factors.single.tissue.sh ${line} &>${RAREVARDIR}/logs/calc.PEER.factors.${line}.log
}

# Calculate PEER factors for each tissue
# Run in parallel requesting 10 cores
export -f factorCalculator
for line in $(cat preprocessing/PEER/gtex_eqtl_tissues.txt); do echo $line; done | parallel --jobs 10 factorCalculator {1}
