#!/usr/bin/env python

import os
import sys
import fileinput
import argparse
import gzip
import numpy as np
import scipy as sp
import random


############## FUNCTIONS

def get_inds(DIR, itype = 'outliers'):
	'''
	Function to retrieve list of individuals with outliers or feature files.
	'''
	inds = np.array(os.listdir(DIR))
	if itype == 'outliers':
		inds = inds[[i for i in range(len(inds)) if 'bed' in inds[i]]]
	elif itype == 'features':
		inds = inds[[i for i in range(len(inds)) if 'GTEX' in inds[i] and 'bed' in inds[i]]]
	else:
		sys.exit('Incorrect itype. Must be features or outliers.')
	inds = [ind[0:9] for ind in inds]
	return inds

def get_inds_from_file(INDS):
	'''
	Function to retrieve list of individuals from a file
	'''
	inds = []
	INDS_fh = open(INDS)
	for line in INDS_fh:
		line = line.rstrip()
		inds.append(line)
	return inds

def make_outlier_dict(OUTLIER_PICKED, individs, threshold, TYPE):
	'''
	Function to make dict of outlier genes.
	Keys will be gene IDs for genes with at least one outlier.
	Values will be dict of outlier IDs for that gene and a ranking metric for each outlier (i.e. FDR). 
	'''
	outlier_dict = {}
	outlier_file = open(OUTLIER_PICKED)
	outlier_file.readline()
	if TYPE == 'Z':
		rank_index = 3
		xins = True
	elif TYPE == 'BF':
		rank_index = 5
		xins = False
	elif TYPE == 'BH':
		rank_index = 6
		xins = False
	else:
		sys.exit('Need to specify rank type. One of BF, BH, or Z.')
	for line in outlier_file:
		line = line.rstrip().split('\t')
		gene = line[0]
		ind = line[1]
		rank = float(line[rank_index])
		if ind in individs:
			if (xins and np.abs(rank) >= threshold) or (not xins and rank <= threshold):
				if gene not in outlier_dict:
					outlier_dict[gene] = {'outliers' : [], 'ranks' : [], 'controls' : []}
				outlier_dict[gene]['outliers'].append(ind)
				outlier_dict[gene]['ranks'].append(str(rank))
	return outlier_dict

def make_test_dict(COUNTS, individs, count_threshold):
	'''
	Function to find list of individuals tested for outliers for each gene.
	'''
	COUNTS_fh = open(COUNTS)
	test_dict = {}
	header = np.array(COUNTS_fh.readline().rstrip().split('\t'))
	indices_to_keep = [0]
	for i in range(1,len(header)):
		if header[i] in individs:
			indices_to_keep.append(i)
	header = header[indices_to_keep[1:]]
	for counts_line in COUNTS_fh:
		counts_line = np.array(counts_line.rstrip().split('\t'))[indices_to_keep]
		gene = counts_line[0]
		counts = np.array([float(x) for x in counts_line[1:]])
		tested = counts >= count_threshold
		test_dict[gene] = header[tested]
	return test_dict

def find_controls(FEATURE_DIR, END, individs, outlier_dict, test_dict):
	'''
	Function to assign controls matched by gene.
	Modifying the outlier dict produced above. Finds all possible controls for a 
	given gene ensuring the controls were tested for outliers.
	Finally, inverts the outlier dict to yield a dict with individual IDs as 
	keys and lists of genes to extract features for for that individual.
	'''
	for ind in individs:
		feature_file = open(FEATURE_DIR + '/' + ind + END)
                feature_file.readline() # skip header
		for line in feature_file:
			line = line.rstrip().split('\t')
			gene = line[0]
			if gene in outlier_dict and ind in test_dict[gene]:
				if ind not in outlier_dict[gene]['outliers']:
					outlier_dict[gene]['controls'].append(ind)
	
	assigned = {}
	for gene in outlier_dict:
		outliers = outlier_dict[gene]['outliers']
		nout = len(outliers)
		controls = outlier_dict[gene]['controls']
		ncon = len(controls)
		statuses = np.concatenate([np.repeat(['1'], nout), np.repeat(['0'], ncon)])
		ranks = np.concatenate([np.array(outlier_dict[gene]['ranks']), np.repeat(outlier_dict[gene]['ranks'][0], ncon)])
		ids = list(outliers) + list(controls)
		for i in range(len(ids)):
			if ids[i] not in assigned:
				assigned[ids[i]] = {'genes' : [], 'statuses' : [], 'ranks' : []}
			assigned[ids[i]]['genes'].append(gene)
			assigned[ids[i]]['statuses'].append(statuses[i])
			assigned[ids[i]]['ranks'].append(ranks[i])
	return assigned

def find_features(assigned, FEATURE_DIR, END, OUT):
	'''
	Function to output features, ranks and outlier statuses for controls and outliers at each gene.
	'''
	OUT_fh = open(OUT, 'w')
	individs = assigned.keys()
	feature_file = open(FEATURE_DIR + '/' + individs[0] + END)
	header = feature_file.readline().rstrip() + '\t' + 'ID' + '\t' + 'RANK' + '\t' + 'Y' + '\n'
	OUT_fh.write(header)
	out_dict = {}
	for ind in individs:
		genes = assigned[ind]['genes']
		statuses = assigned[ind]['statuses']
		ranks = assigned[ind]['ranks']
		feature_file = open(FEATURE_DIR + '/' + ind + END)
		for line in feature_file:
			line = line.rstrip()
			split_line = line.split('\t')
			gene = split_line[0]
			if gene in genes:
				index = genes.index(gene)
				rank = ranks[index]
				status = statuses[index]
				out_line = line + '\t' + ind + '\t' + rank + '\t' + status + '\n'
				if gene not in out_dict:
					out_dict[gene] = {'outliers' : [], 'controls' : []}
				if status == '1':
					out_dict[gene]['outliers'].append(out_line)
				else:
					out_dict[gene]['controls'].append(out_line)
	for gene in out_dict:
		nout = len(out_dict[gene]['outliers'])
		ncon = len(out_dict[gene]['controls'])
		if nout >= 1 and ncon >= 1:
			feature_lines = out_dict[gene]['outliers'] + out_dict[gene]['controls']
			for line in feature_lines:
				OUT_fh.write(line)
	OUT_fh.close()

if __name__ == "__main__":
	
	############## MAIN
	USAGE = """
	Picks outliers and matched controls from full feature sets compiled for each individual.
    """

	parser = argparse.ArgumentParser(description = USAGE)
	parser.add_argument('--OUTLIER_PICKED', dest = 'OUTLIER_PICKED', required = True, help = 'outlier file')
	parser.add_argument('--COUNTS', dest = 'COUNTS', required = True, help = 'counts file')
	parser.add_argument('--FEATURE_DIR', dest = 'FEATURE_DIR', required = True, help = 'directory holding features files')
	parser.add_argument('--OUT', dest = 'OUT', required = True, help = 'output file')
	parser.add_argument('--threshold', dest = 'threshold', default = 0.01, help = 'threshold for calling significant outliers')
	parser.add_argument('--count_threshold', dest = 'count_threshold', default = 5, help = 'threshold for number of observed tissues')
	parser.add_argument('--type', dest = 'TYPE', default = 'BH', help = 'rank type one of either Z, BH, or BF.')
	parser.add_argument('--INDS', dest = 'INDS', required = True, help = 'list of individuals to consider for feature construction')
	parser.add_argument('--END', dest = 'END', required = True, help = 'feature file suffix')
	
	args = parser.parse_args()

	INDS = args.INDS
	OUTLIER_PICKED = args.OUTLIER_PICKED
	FEATURE_DIR = args.FEATURE_DIR
	OUT = args.OUT 
	threshold = float(args.threshold)
	TYPE = args.TYPE
	END = args.END
	count_threshold = float(args.count_threshold)
	COUNTS = args.COUNTS

	## Get list of individuals to consider from the gene level feature files
	individs = get_inds_from_file(INDS)

	## Make dict of individuals tested for outliers
	test_dict = make_test_dict(COUNTS, individs, count_threshold)

	## Get list of outlier individuals
	outlier_dict = make_outlier_dict(OUTLIER_PICKED, individs, threshold, TYPE)

	## Pick matched controls
	assigned = find_controls(FEATURE_DIR, END, individs, outlier_dict, test_dict)

	## Output features for outliers and controls
	find_features(assigned, FEATURE_DIR, END, OUT)
