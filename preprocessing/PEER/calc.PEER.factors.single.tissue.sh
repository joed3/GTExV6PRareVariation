#!/bin/bash

set -o nounset -o errexit -o pipefail

PEER_DIR=${RAREVARDIR}/preprocessing/PEER
TraitsFileName=${PEER_DIR}/${1}.rpkm.log2.ztrans.txt
CovsFileName=${RAREVARDIR}/data/eqtl_data/eQTLInputFiles/covariates/${1}_Analysis.covariates.txt
MaxFactorsN=$(grep InferredCov ${CovsFileName} | wc -l)
MaxIterations=10000
BoundTol=0.001
VarTol=0.00001
e_pa=0.1
e_pb=10.
a_pa=0.001
a_pb=0.1
OutDir=${PEER_DIR}/${1}_Factors"$MaxFactorsN"
tissue=${1}

echo "R -f calc.PEER.factors.core.R --vanilla --slave --args $TraitsFileName $MaxFactorsN $MaxIterations $BoundTol $VarTol $e_pa $e_pb $a_pa $a_pb $OutDir $tissue"
R -f ${PEER_DIR}/calc.PEER.factors.core.R --vanilla --slave --args $TraitsFileName $MaxFactorsN $MaxIterations $BoundTol $VarTol $e_pa $e_pb $a_pa $a_pb $OutDir $tissue

echo "calc.PEER.factors.single.tissue.sh DONE"


