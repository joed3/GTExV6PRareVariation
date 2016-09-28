#!/bin/bash

set -o nounset -o errexit -o pipefail

## run add_features_variant_beds.sh in batches

nproc=15
dir=${RAREVARDIR}/features

parallel --jobs $nproc --xapply ./add_features_variant_beds.sh ::: ${dir}/variantBeds/individuals/withCADD/*_SNPs.CADD.bed.gz ::: ${dir}/variantBeds/1KG/SNPs.1kg.AF.bed.gz ::: ${dir}/bySite

parallel --jobs $nproc --xapply ./add_features_variant_beds.sh ::: ${dir}/variantBeds/individuals/*_indels.bed.gz ::: ${dir}/variantBeds/1KG/indels.1kg.AF.bed.gz ::: ${dir}/bySite

wait
