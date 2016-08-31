#!/usr/bin/env python

## Reads in the raw Pouya Motifs bed file and makes some summary statistics

import sys
import numpy as np
import gzip

infile = sys.argv[1]
outfile = sys.argv[2]

## Read in Pouya motifs dataset and make dict with keys as factors and motifs as values
	## sub-dict will have motifs as keys and lengths and number of occurrences as values
tf_dict = {}
POUYA = gzip.open(infile, 'rb')

for line in POUYA:
	chrom, start, stop, motif, strand, factor = line.rstrip().split('\t')
	start = float(start)
	stop = float(stop)
	length = stop - start
	if factor not in tf_dict:
		tf_dict[factor] = {}
	if motif not in tf_dict[factor]:
		tf_dict[factor][motif] = {'avg_length' : 0, 'counts' : 0}
	tf_dict[factor][motif]['avg_length'] += length
	tf_dict[factor][motif]['counts'] += 1

POUYA.close()

## Calculate average length for each motif 
## Write output to file for easy plotting with R
OUT = open(outfile, 'w')
OUT.write('TF\tMOTIF\tLENGTH\tCOUNTS\n')
for factor in tf_dict:
	for motif in tf_dict[factor]:
		tf_dict[factor][motif]['avg_length'] = float(tf_dict[factor][motif]['avg_length']) / tf_dict[factor][motif]['counts']
		OUT.write('\t'.join([factor, motif, str(tf_dict[factor][motif]['avg_length']), str(tf_dict[factor][motif]['counts'])]) + '\n')
OUT.close()
