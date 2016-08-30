#!/usr/bin/env python

# takes in gtf file as input and extracts protein-coding and lincRNA exons
# while doing so, pads internal exons by 5 bp to capture the canonical splice sites
# outputs a bedfile to stdout

from __future__ import print_function
import pybedtools
import argparse
import sys

class Transcript(object):
    def __init__(self, chrom, start, end):
        self._chr = chrom
        self._starts = [start]
        self._ends = [end]
    
    def __repr__(self):
        value = "Chrom " + self._chr + "\n"
        value = value + "Starts " + repr(self._starts) + "\n"
        value = value + "Ends " + repr(self._ends) + "\n"
        return value

    def printPadded(self, pad = 5):
        """
        Print out transcript, one line for each exon, in bed format.
        Pad each internal exon with the given number of base pairs on each side.
        So pad the extremities of every exon except the the start of the first exon 
        and the end of the last.
        Outputs to stdout. Returns nothing.
        """
        # if there is a single exon: just print it
        if len(self._starts) == 1:
            print(elements2bed(self._chr, self._starts[0], self._ends[0]))
            return
        # first sort start and end coordinates
        # assuming that the transcripts are well defined in that none of the exons overlap
        self._starts.sort()
        self._ends.sort()
        pairs = zip(self._starts, self._ends)
        # first deal with last and first exons
        last = pairs.pop()
        first = pairs.pop(0)
        print(elements2bed(self._chr, first[0], (first[1] + pad)))
        print(elements2bed(self._chr, (last[0] - pad), last[1]))
        for p in pairs:
            print(elements2bed(self._chr, (p[0] - pad), (p[1] + pad)))

    def addExon(self, start, end):
        self._starts.append(start)
        self._ends.append(end)

def elements2bed(chrom, start, stop):
    return '\t'.join([chrom, str(start), str(stop)])

def main():
    # process arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help = 'Path to gtf file to process.')
    args = parser.parse_args()
    
    # read in gtf file as BedTool
    gtf = pybedtools.BedTool(args.gtf)
    
    # filter gtf to only contain exons from pretein-coding and lincRNA genes
    filtered = gtf.filter(lambda i: i[2] == 'exon' and i.attrs['gene_type'] in ['protein_coding','lincRNA'])

    # make no assumptions about how lines are grouped together
    # process the filtered gtf lines into a dict {transcript: transcript object}
    transDict = dict()
    for line in filtered:
        transcript = line.attrs['transcript_id']
        if transcript not in transDict:
            transDict[transcript] = Transcript(line.chrom, line.start, line.stop)
        else:
            transDict[transcript].addExon(line.start, line.stop)

    # Print out padded exons in bed format
    for t in transDict.values():
        t.printPadded()
            

if __name__ == '__main__':
    main()
