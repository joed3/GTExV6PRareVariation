#!/usr/bin/env python

from __future__ import division
import sys
import argparse
from collections import defaultdict
import pandas as pd
import numpy as np
import os

###
### Command line argument parser
###
parser = argparse.ArgumentParser(description="Build features file/response file for project.")
parser.add_argument('--features', action='store', default="features.txt", type=str, help='the file to write feature rows to')

### function to remove NAs
def rmNA(values):
        filtered = []
        for val in values:
                if val !="NA":
                        filtered.append(val)
        return filtered

### function to return max of given column (returns NA if empty or no observed values)
def get_max(x, colname, data_type=float):
        vals = rmNA(pd.Series(x[colname]))
        if len(vals) == 0: return 'NA'
        maxval = max(map(data_type, vals))
        if data_type==int: return maxval
        return round(maxval,3)

### function to see whether any of the values are non-zero (returns NA if no variants, but returns 0 if not observed)
def get_any(x, colname):
        vals = rmNA(pd.Series(x[colname]))
        if len(vals) == 0: return 'NA'
        return 1 if np.any(map(lambda a: int(a) > 0, vals)) else 0

###
### Define the functions that calculate individual features
###
def gene_id(x):
	return x['geneID'][0]

def n_variants(x):
	return x.shape[0]

def any_tfbs(x):
	return get_any(x, 'nTFBS')

def any_enhancer(x):
        return get_any(x, 'nenh')

def any_promoter(x):
        return get_any(x, 'nprom')

def any_dyadic(x):
        return get_any(x, 'ndyadic')

def any_tfbs_CADD(x):
	return get_any(x, 'nTFBSCADD')

def any_tfbs_CADD_peaks(x):
	return get_any(x, 'TFBSCADDPeaks')

def max_tfbs(x):
        return get_max(x, 'nTFBS', data_type=int)

def max_enhancer(x):
        return get_max(x, 'nenh', data_type=int)

def max_promoter(x):
        return get_max(x, 'nprom', data_type=int)

def max_dyadic(x):
        return get_max(x, 'ndyadic', data_type=int)

def max_tfbs_CADD(x):
        return get_max(x, 'nTFBSCADD', data_type=int)

def max_tfbs_CADD_peaks(x):
        return get_max(x, 'TFBSCADDPeaks', data_type=int)

def max_tfbs_CADD_peaks_max(x):
        return get_max(x, 'TFBSCADDPeaksMax')

def max_cadd(x):
        return get_max(x, 'CADD')

def max_pri_phylop(x):
	return get_max(x, 'priPhyloP')

def max_mam_phylop(x):
	return get_max(x, 'mamPhyloP')

def max_ver_phylop(x):
	return get_max(x, 'verPhyloP')

def max_pri_phastCons(x):
	return get_max(x, 'priPhCons')

def max_mam_phastCons(x):
	return get_max(x, 'mamPhCons')

def max_ver_phastCons(x):
	return get_max(x, 'verPhCons')

def max_fitCons(x):
	return get_max(x, 'fitCons')

def max_GerpN(x):
	return get_max(x, 'GerpN')

def max_GerpS(x):
	return get_max(x, 'GerpS')

def max_GerpRS(x):
	return get_max(x, 'GerpRS')

def max_CpG(x):
        return get_max(x, 'CpG')

def closest_variant(x):
        n = n_variants(x)
        if n==0: return 'NA'
        tsspos = int(x['TSS1'][0])
	positions = np.array(map(int, np.array(pd.Series(x['Pos1']))))
	return np.amin(np.abs(positions - tsspos))

###
### Order the features as we want them in the output file
###
features = []
features.append(gene_id)
features.append(n_variants)
features.append(any_tfbs)
features.append(any_enhancer)
features.append(any_promoter)
features.append(any_dyadic)
features.append(any_tfbs_CADD)
features.append(any_tfbs_CADD_peaks)
features.append(max_tfbs)
features.append(max_enhancer)
features.append(max_promoter)
features.append(max_dyadic)
features.append(max_tfbs_CADD)
features.append(max_tfbs_CADD_peaks)
features.append(max_tfbs_CADD_peaks_max)
features.append(max_cadd)
features.append(max_pri_phylop)
features.append(max_mam_phylop)
features.append(max_ver_phylop)
features.append(max_pri_phastCons)
features.append(max_mam_phastCons)
features.append(max_ver_phastCons)
features.append(max_fitCons)
features.append(max_GerpN)
features.append(max_GerpS)
features.append(max_GerpRS)
features.append(max_CpG)
features.append(closest_variant)

###
### Configure input
###
columnNames = ['chromosome','TSS0','TSS1','geneID','variantChromosome','Pos0','Pos1','MAF','genotype','CADD','CADDnorm','nTFBS','nprom','nenh','ndyadic','CpG','priPhCons','mamPhCons','verPhCons','priPhyloP','mamPhyloP','verPhyloP','GerpN','GerpS','GerpRS','GerpRSpval','fitCons','nTFBSCADD','TFBSCADDPeaks','TFBSCADDPeaksMax','EAS_AF','AMR_AF','AFR_AF','EUR_AF','SAS_AF']

###
### Main function
###
if __name__ == "__main__":

	### Print command line
	print >> sys.stderr, "[Command]", " ".join(sys.argv)

	### Parse command line arguments
	args = parser.parse_args()

        ### open output file and add header
        featuresFile = open(args.features, 'w')
        featuresFile.write("\t".join([f.__name__ for f in features]) + "\n")

	### Store all the lines in appropriate subgroups per gene
	variants_per_gene = defaultdict(list)

	### Process stdin
	for line in sys.stdin:

		### Breakdown line
		fields = line.strip().split('\t')
                
                ### Use gene ID as key
                key = fields[3]

                ### Store variant data with the gene ID
                variants_per_gene[key].append(fields)
	
	### Iterate over collection of genes
	for gene, variants in variants_per_gene.iteritems():

		### Create a dataframe that we can use logically
		data = pd.DataFrame(variants, columns=columnNames)

		### Collapse variants to create features
		collapsed_features = [f(data) for f in features]

		featuresFile.write("\t".join([str(s) for s in collapsed_features]) + '\n')

	### Close output file
	featuresFile.close()
