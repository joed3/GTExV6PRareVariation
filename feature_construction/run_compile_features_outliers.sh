#!/bin/bash

# run get v6 features for median Z
bash compile_features_outliers.sh

# get v6 features for median Z without a threshold on the Z-score (only for rarest MAF bin)
bash compile_features_outliers_nothresh.sh

# get v6 features for singlez (paths are hardcoded in the script)
bash compile_features_outliers_singletissue.sh
