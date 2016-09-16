#!/usr/bin/env python

# Python script to perform multiple testing correction (BF on the gene level) on GTEx Metasoft results

import os
import sys
import fileinput
import argparse
import gzip
import numpy as np
import scipy as sp

#------------- FUNCTIONS

def bfCorrect(META_fh, TISS_fh, OUT_fh):
	META = open(META_fh)
	TISS = open(TISS_fh)
	pnames = []
	mnames = []
	for line in TISS:
		line = line.rstrip().split('\t')
		pnames.append('P.' + line[0])
		mnames.append('M.' + line[0])
	TISS.close()
	OUT = open(OUT_fh, 'w')
	header = META.readline().rstrip().split('\t')
	header = ['SNP'] + header
	header[1] = 'GENE'
	header[2] = 'N.TISSUES'
	header.pop()
	header.pop()
	header = header + pnames + mnames + ['pvalRE2.BF', 'TESTS']
	OUT.write('\t'.join(header) + '\n')
	bfDict = {}
	for line in META:
		line = line.rstrip().split('\t')
		snp, gene = line[0].split(',')
		if line[8] != 'NA':
			pvalRE2 = float(line[8])
			if gene not in bfDict:
				bfDict[gene] = {'snp' : snp, 'pvalRE2' : pvalRE2, 'tests' : 0, 'line' : '\t'.join([snp, gene] + line[1:])}
			if pvalRE2 < bfDict[gene]['pvalRE2']:
				bfDict[gene]['snp'] = snp
				bfDict[gene]['pvalRE2'] = pvalRE2
				bfDict[gene]['line'] = '\t'.join([snp, gene] + line[1:])
			bfDict[gene]['tests'] += 1
	META.close()
	for gene in bfDict:
		bf = np.min([bfDict[gene]['pvalRE2'] * bfDict[gene]['tests'], 1])
		OUT.write(bfDict[gene]['line'] + '\t' + str(bf) + '\t' + str(bfDict[gene]['tests']) + '\n')
	OUT.close()

#------------- MAIN

USAGE = """
	Performs BF correction on the gene level on Metasoft results from GTEx.
    """

parser = argparse.ArgumentParser(description = USAGE)
parser.add_argument('--META', dest = 'META', required = True, help = 'Metasoft results from GTEx')
parser.add_argument('--TISS', dest = 'TISS', required = True, help = 'Metasoft tissue order')
parser.add_argument('--OUT', dest = 'OUT', required = True, help = 'output file to hold BF corrected data')

args = parser.parse_args()

META_fh = args.META
TISS_fh = args.TISS
OUT_fh = args.OUT

# Perform BF correction 
bfCorrect(META_fh, TISS_fh, OUT_fh)
