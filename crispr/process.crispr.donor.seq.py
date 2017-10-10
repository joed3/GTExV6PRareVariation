import os
import sys
import fileinput
import argparse
import gzip
import numpy as np
import scipy as sp
import re

#----------- FUNCTIONS

# Function to process FASTA file
def processFasta(FASTA):
	fasta = {}
	FASTA = open(FASTA)
	ID = ''
	for line in FASTA:
		line = line.rstrip()
		if line[0] == '>':
			ID = line
			fasta[line] = {'seq' : '', 'ref' : '', 'alt' : ''}
		else:
			fasta[ID]['seq'] += line
	FASTA.close()
	return fasta

# Function to process the VCF file
def processVCF(VCF, fasta, radius):
	VCF = open(VCF)
	VCF.readline()
	idList = []
	for line in VCF:
		line = line.rstrip().split('\t')
		chrom = line[0]
		pos = int(line[1])
		ref = line[3]
		alt = line[4]
		ID = '>chr' + chrom + ':' + str(pos - radius) + '-' + str(pos + radius)
		fasta[ID]['ref'] = ref
		fasta[ID]['alt'] = alt
		idList.append(ID)
	VCF.close()
	return idList

# Generate donor sequence
def generateDonorSeq(fasta, radius, idList, OUT):
	OUT = open(OUT, 'w')
	for ID in idList:
		refSeq = fasta[ID]['seq']
		altSeq = fasta[ID]['seq'][:radius] + fasta[ID]['alt'] + fasta[ID]['seq'][radius + 1:] 
		# FASTA format has at most 80 characters per line
		n = 80
		refSeq = [refSeq[i:i+n] for i in range(0, len(refSeq), n)]		
		OUT.write(ID + ':' + fasta[ID]['ref'] + '\n')
		for i in refSeq:
			OUT.write(i.upper() + '\n')
		
		altSeq = [altSeq[i:i+n] for i in range(0, len(altSeq), n)]
		OUT.write(ID + ':' + fasta[ID]['alt'] + '\n')
		for i in altSeq:
			OUT.write(i.upper() + '\n')
	OUT.close()

#----------- MAIN
if __name__ == "__main__":
	
	USAGE = """
	Generate donor sequences for CRISPR. One sequence with the ref and one with the alt for each variant. FASTA format.
    """

	parser = argparse.ArgumentParser(description = USAGE)
	parser.add_argument('--VCF', dest = 'VCF', required = True, help = 'input VCF')
	parser.add_argument('--FASTA', dest = 'FASTA', required = True, help = 'input FASTA file with raw sequence to process')
	parser.add_argument('--OUT', dest = 'OUT', required = True, help = 'output file')
	parser.add_argument('--radius', dest = 'radius', required = False, help = 'length of sequence to the left and right of the variant', default = 49)
	
	args = parser.parse_args()

	VCF = args.VCF
	FASTA = args.FASTA
	OUT = args.OUT
	radius = int(args.radius)

	# Process the raw donor sequence FASTA
	fasta = processFasta(FASTA)

	# Process the VCF
	idList = processVCF(VCF, fasta, radius)

	# Generate the donor sequences
	generateDonorSeq(fasta, radius, idList, OUT)

