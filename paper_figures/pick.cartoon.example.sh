#!/bin/bash

dir=$RAREVARDIR

cat ${dir}/data/outliers_medz_picked.txt | awk '$3==5' > possible.txt

# subset the normalized expression file to these example genes
expr=${dir}/preprocessing/gtex_2015-01-12_normalized_expression.txt
head -n1 $expr > expr.subset.by.genes.txt
cat possible.txt | cut -f1 | grep -f - $expr >> expr.subset.by.genes.txt


