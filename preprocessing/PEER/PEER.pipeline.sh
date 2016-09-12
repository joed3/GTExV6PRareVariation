#!/bin/bash

# PEER correction pipeline for outlier detection and rare variant analysis

PEER_DIR=${1}

# Make RPKM and read count matrices for each tissue
R CMD BATCH --no-save ${PEER_DIR}/preprocess.expr.data.R &>${RAREVARDIR}/logs/preprocess.expr.data.log

# Estimate PEER factors
bash ${PEER_DIR}/calc.PEER.factors.all.tissues.sh ${PEER_DIR}

# Remove residuals
bash ${PEER_DIR}/calc.residuals.sh
