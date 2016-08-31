#!/usr/bin/env python

# author: Emily Tsang

# Takes allele frequencies for all populations extracted from VCF (from the command line)
# sums allele frequencies for different alleles at the same position

# outputs in bed format to stdout

import sys

prevChrom = ''
prevPos = ''
prevFreqs = [0,0,0,0,0]

def maf(altAF):
    refAF = 1 - altAF
    return min(altAF, refAF)

for line in sys.stdin:
    line = line.strip().split()
    # skip header line
    if line[0] == "CHROM":
        continue

    chrom = line[0]
    pos = line[1]
    freqs = line[4:]

    # if multiple frequencies per column (because of multiple alleles), add them
    freqs = [reduce(lambda a,b: a+b, map(float,f.split(','))) for f in freqs]
    freqs = [maf(f) for f in freqs]

    # add up frequencies as long as it's the same position
    if chrom == prevChrom and pos == prevPos:
        prevFreqs = [max(a,b) for a,b in zip(prevFreqs,freqs)]
    # if new position output previous position
    else:
        if prevChrom!='':
            sys.stdout.write("\t".join(["chr"+prevChrom, str(int(prevPos)-1), prevPos] + [str(round(f,4)) for f in prevFreqs]) + '\n')
        prevFreqs = freqs
        prevChrom = chrom
        prevPos = pos

# print last line
sys.stdout.write("\t".join(["chr"+prevChrom, str(int(prevPos)-1), prevPos] + [str(round(f,4)) for f in prevFreqs]) + '\n')
