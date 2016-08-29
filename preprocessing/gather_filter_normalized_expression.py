#!/usr/bin/env python

# Takes PEER-corrected expression data 
# and combines it into a flat file.
# The output file is formatted as follows:
# 
#	tissue gene ind1_expr ind2_expr ...
#
# Missing values are coded as NAs
# Gets tissues and individual IDs from file

import numpy as np
from operator import itemgetter

dir = '../preprocessing/'
tissueNamesFile = dir + 'gtex_2015-01-12_tissues_all_normalized_samples.txt'
individualsFile = dir + 'gtex_2015-01-12_individuals_all_normalized_samples.txt'
outfile = dir + 'gtex_2015-01-12_normalized_expression.txt'

exprdir = dir + 'PEER/'

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

# prepare output file
out = open(outfile, 'w')
out.write('\t'.join(['Tissue','Gene'] + individuals) + '\n')

# process tissues one at a time 
for tissue in tissues:
	filename = exprdir + tissue + '.peer.ztrans.txt'

	# read in header and figure out which columns to keep
	expr = open(filename, 'r')
	headerList = expr.readline().strip().split()
	subjectIDs = headerList[1:]
	cols2keep = [(i+1,subject) for (i,subject) in enumerate(subjectIDs) if subject in individuals]
        # sort by subject ID (because brain_cerebellum wasn't sorted for some reason)
        sorted2keep = sorted(cols2keep, key=itemgetter(1))
	cols2keep = [0] + [s2k[0] for s2k in sorted2keep]
	index2add = [1] + [i+2 for (i,subject) in enumerate(individuals) if subject in subjectIDs]
	ngenes = sum(1 for _ in expr)
	expr.close()

	# make sure that individual IDs line up
	orig = [(['Id'] + individuals)[i-1] for i in index2add ]
	new = [headerList[i] for i in cols2keep]
	assert len(new) == sum([o==n for o,n in zip(orig,new)])

	# initialize numpy array to contain final data
	# fill with NAs
	padded = np.full((ngenes,len(individuals)+2), 'NA', dtype='|S40')
	# add column with tissue
	padded[:,0] = tissue

	# read data into numpy array
	values = np.loadtxt(filename, dtype=str, skiprows=1)
	# remove unwanted columns
	values = values[:, cols2keep]
	# put the values into the right slots in the preallocated array
	padded[:,index2add] = values
	# print matrix to file
	np.savetxt(out, padded, fmt='%s', delimiter='\t', newline='\n')

out.close()
