#!/bin/bash

# run get v6 features for median Z with top 0 PEER factors removed
bash feature_construction/compile_features_outliers_peerless.sh .peer.top0.txt &> ${RAREVARDIR}/logs/run_compile_features_outliers_peerless.top0.log &

# run get v6 features for median Z with top 5 PEER factors removed
bash feature_construction/compile_features_outliers_peerless.sh .peer.top5.txt &> ${RAREVARDIR}/logs/run_compile_features_outliers_peerless.top5log &
