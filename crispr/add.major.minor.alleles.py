import os
import sys
import fileinput
import argparse
import gzip
import numpy as np
import scipy as sp
import re

#----------- FUNCTIONS


#----------- MAIN
if __name__ == "__main__":
	
	USAGE = """
	Add info about major and minor alleles to annotated CRISPR variants.
    """

	parser = argparse.ArgumentParser(description = USAGE)
	parser.add_argument('--IN', dest = 'IN', required = True, help = 'input VCF')
	parser.add_argument('--OUT', dest = 'OUT', required = True, help = 'output VCF')
	parser.add_argument('--FRQ', dest = 'FRQ', required = True, help = 'frequency data')

	args = parser.parse_args()

	IN = args.IN
	OUT = args.OUT
	FRQ = args.FRQ

	IN = open(IN)
	OUT = open(OUT, 'w')
	FRQ = open(FRQ)

	alleles = {}
	FRQ.readline()
	for line in FRQ:
		line = line.rstrip().split('\t')
		chrom = line[0]
		pos = line[1]
		if chrom not in alleles:
			alleles[chrom] = {}
		alleles[chrom][pos] = {}
		for allele in line[4:]:
			nuc, frq = allele.split(':')
			alleles[chrom][pos][nuc] = float(frq)
	FRQ.close()

	header = IN.readline().rstrip()
	OUT.write(header + '\t' + 'Major' + '\n')
	for line in IN:
		line = line.rstrip().split('\t')
		chrom = line[0]
		pos = line[1]
		ref = line[3]
		alt = line[4]
		if alleles[chrom][pos][ref] >= 0.5:
			major = ref
		else:
			major = alt
		OUT.write('\t'.join(line + [major]) + '\n')



