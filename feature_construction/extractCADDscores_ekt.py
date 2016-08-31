#!/usr/bin/env python

"""
Modeled after script provided on CADD website written by Martin Kircher

Reads a bed file from std in with the following columns:
chr pos0 pos1 MAF genotype allele1 allele2

Gets the CADD score for each allele (only exists the non reference ones) and records the max of the two

Prints:
chr pos0 pos1 MAF genotype CADD
CADD will be NA if the individual is homozygous for the reference

Run as follows:
cat input.bed | python extractCADDscores_ekt.py 
"""

import pysam
import sys
import os

caddDir = os.environ["CADD_DIR"]
filename = caddDir + "/mnt/lab_data/montgomery/shared/CADD/whole_genome_SNVs.tsv.gz"

# BED FIELDS
fchr = 0
fpos1 = 2
fallele1 = 5
fallele2 = 6

# CADD FIELDS
alt = 3
rawscore = 4
phredscore = 5

sys.stderr.write("Opening %s...\n"%(filename))
regionTabix = pysam.Tabixfile(filename,'r')

for line in sys.stdin:
  fields = line.rstrip().split('\t')
  chrom = fields[fchr]
  # strip "chr"
  chrom = chrom[3:]
  pos = int(fields[fpos1])
  alleles = set([fields[fallele1],fields[fallele2]])

  raw = []
  phred = []
  for allele in alleles:
    for CADDline in regionTabix.fetch(chrom,pos-1,pos):
      caddFields = CADDline.rstrip().split('\t')
      if (caddFields[alt] == allele):
        raw.append(float(caddFields[rawscore]))
        phred.append(float(caddFields[phredscore]))
        break
  
  # print max of scores or NA if scores is empty
  if not raw:
    sys.stdout.write("\t".join(fields[0:5] + ['NA','NA']) + '\n')
  else:
    maxRaw = max(raw) # need max in case of multiple alternate alleles
    maxPhred = max(phred)
    sys.stdout.write("\t".join(fields[0:5] + [str(maxRaw),str(maxPhred)]) + '\n')
