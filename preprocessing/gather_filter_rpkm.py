#!/usr/bin/env python

# Make flat file with RPKM
# With columns: 
#	Tissue Gene Ind1 Ind2 ...
# Missing values are coded as NAs

import os
import sys
import fileinput
import argparse
import gzip
import numpy as np
import scipy as sp

#------------ FUNCTIONS

def getInputOutputIndices(indivs, columnDict):
	'''function that gets input and output column indices
	args: list of individuals to include in output (ordered),
		  dict individual as key and input column index as value
	return: a tuple of two lists; 
			the first is the columns to grab from the input matrix
			the second is the columns to insert them into in the output matrix
	'''
	outCols = [1] # for gene id
	inCols = [0] # column with ensg
	# iterate through individuals and establish if they have data for that tissue
	for index, indiv in enumerate(indivs):
		if indiv in columnDict:
			outCols.append(index+2) # +2 because of the tissue and the gene columns
			inCols.append(columnDict[indiv])
	return (inCols, outCols)

def getIndiv(sample):
        splitsamp = sample.split('-')
        return splitsamp[0] + '-' + splitsamp[1]


#----------- MAIN

iodir = os.environ["RAREVARDIR"]
iodir = iodir + '/preprocessing/'
tissueNamesFile = iodir + 'gtex_2015-01-12_tissues_all_normalized_samples.txt'
individualsFile = iodir + 'gtex_2015-01-12_individuals_all_normalized_samples.txt'
rpkmDir = iodir + 'PEER/'
outfile = iodir + 'gtex_2015-01-12_rpkm.txt'

# read IDs into a list
# add 'GTEX-' as a prefix and sort by id
ind = open(individualsFile, 'r')
individuals = [i.strip() for i in ind.readlines()]
individuals.sort()
ind.close()

# read in tissues to process in a list
# sort it
tis = open(tissueNamesFile, 'r')
tissues = [t.strip() for t in tis.readlines()]
tissues.sort()
tis.close()

# read in RPKM file for each tissue
# reorder output and write out with padded values for missing individuals
out = open(outfile, 'w')
out.write('\t'.join(['Tissue', 'Gene'] + individuals) + '\n')
out.close()

for tissue in tissues:
	print tissue
	rpkmFile = rpkmDir + tissue + '.rpkm.txt'
	values = np.loadtxt(rpkmFile, dtype=str)
	ngenes = values.shape[0] - 1
	header = values[0, 1:]
	padded = np.full((ngenes,len(individuals)+2), 'NA', dtype='|S40')
	padded[:, 0] = tissue
	padded[:, 1] = values[1:, 0]
	for i in range(len(header)):
		if header[i] in individuals:
			outIndex = individuals.index(header[i]) + 2
			padded[:, outIndex] = values[1:, i + 1]
	out = open(outfile, 'a')
	np.savetxt(out, padded, fmt='%s', delimiter='\t', newline='\n')
	out.close()
