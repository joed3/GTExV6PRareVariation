#!/bin/bash

set -o nounset -o errexit -o pipefail

# set number of processes
nproc=25

# get CADD scores for each of the individual SNP bed files

if [ $# -ne 1 ]; then
    echo "usage: compileCADDscores.sh <bed file directory>"
    exit
fi

indir=${1}/individuals # compatible with outdir of vcf2bedfiles_helper_processVCFtoolsOutput.sh
outdir=${indir}/withCADD

# make output directory if it doesn't exist
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

export indir

date
echo

parallel --jobs $nproc "echo {}; cat {} | ./extractCADDscores_ekt.py > {.}.CADD.bed" ::: ${indir}/*_SNPs.bed

# mv all created files to withCADD directory
mv ${indir}/*CADD.bed $outdir

echo
date

echo done
