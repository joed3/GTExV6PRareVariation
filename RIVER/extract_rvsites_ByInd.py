#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
: Author - Yungil Kim 
: E-mail - ipw012@gmail.com

: Description

This script reads targeted regions per gene (chrom, start, stop, gencode_idx, gene_name, TSS). For each gene, search the candidate sites having MAF < 0.01 for GTEx data
and EUR population from 1k genome by navigating AF vcf files. Then, for each inidividual, write a text file having targeted sites (chrom, pos, ref, alt, gencode_idx, gene_name, major==ref?, # var, DistTSS)

cd
source ./.bashrc
cat ${RAREVARDIR}/RIVER/data/rvsite/region.tss10k.txt | ${RAREVARDIR}/RIVER/extract_rvsites_ByInd.py -n 1 --id $ID --WGSvcf_in ${RAREVARDIR}/data/wgs/a_filtered_and_compressed_GTEx_WGS_vcf_file --GTExvcf_in ${RAREVARDIR}/data/wgs/a_compressed_GTEx_allele_frequency_file --EURvcf_in ${RAREVARDIR}/data/wgs/1KG/a_compressed_EUR_allele_freq_vcf_file --site_out ${RAREVARDIR}/RIVER/data/score/indiv/${ID}.txt

"""

import sys, os, re
import pysam
import numpy as np
from optparse import OptionParser
import time
start_time = time.time()

parser = OptionParser()
parser.add_option("-n", type="int", dest="num")
parser.add_option("--id",action="store", type="string", dest="ind_id")
parser.add_option("--WGSvcf_in", dest="WGSvcf_in")    # WGS data
parser.add_option("--GTExvcf_in", dest="GTExvcf_in")    # GTEx allele frequencies
parser.add_option("--EURvcf_in", dest="EURvcf_in")
parser.add_option("--site_out", dest="site_out")  # chrom, pos, ref, alt, gencode_idx, gene_name, major==ref?, # var, DistTSS

(options, args) = parser.parse_args()

## uploading input files
if os.path.exists(options.WGSvcf_in) and os.path.exists(options.WGSvcf_in+".tbi"):
    WGSTabix = pysam.Tabixfile(options.WGSvcf_in,'r') # WGS.vcf
if os.path.exists(options.GTExvcf_in) and os.path.exists(options.GTExvcf_in+".tbi"):
    GTExTabix = pysam.Tabixfile(options.GTExvcf_in,'r') # GTEx.vcf
if os.path.exists(options.EURvcf_in) and os.path.exists(options.EURvcf_in+".tbi"):
    EURTabix = pysam.Tabixfile(options.EURvcf_in,'r') # EUR.vcf


WGS_header = []
with open(options.WGSvcf_in,"r") as f:
    for line in f:
        eachcol = line.rstrip().split('\t')
        if eachcol[0].startswith("#CHROM") == 1:
            for colnames in eachcol:
                WGS_header.append(colnames)
        break

dic_WGS = {}
for indivs in WGS_header:
    dic_WGS[indivs] = WGS_header.index(indivs)

site_out = open(options.site_out,'w')

for target_region in sys.stdin: # targeted regions per gene
    region_info = target_region.rstrip().split('\t')
    gencode_idx = region_info[3]
    gene_name = region_info[4]
    tss_pos = region_info[5]

    for poss_var in WGSTabix.fetch(region_info[0],int(region_info[1]),int(region_info[2])):
        fields_WGSvcf = poss_var.rstrip().split('\t')
        chrom = fields_WGSvcf[0]
        pos = int(fields_WGSvcf[1])

        gt_ind = fields_WGSvcf[int(dic_WGS[options.ind_id])]
        gt_allele = re.split('[/:]',gt_ind) # individual-specific indices of alleles based on ref and alt

        for target_site in GTExTabix.fetch(chrom,pos-1,pos): # AFs in GTEx European subjects
            fields_gtex = target_site.rstrip().split('\t')
            
            count_allele = 0; rv_allele = []; dic_allele_idx = {}; dic_allele_af = {}
            for gtex_allele in fields_gtex[4:]:
                temp_allele = gtex_allele.split(':')
                dic_allele_idx[temp_allele[0]] = count_allele; count_allele += 1
                dic_allele_af[temp_allele[0]] = temp_allele[1]
                if (float(temp_allele[1]) < 0.01) and (float(temp_allele[1]) > 0): # GTEx rv
                    rv_allele.append(temp_allele[0]);

        if len(rv_allele) > 0:
            major_allele = dic_allele_idx.keys()[dic_allele_af.values().index(sorted(dic_allele_af.values())[-1])]

            set_allele = [fields_WGSvcf[3]]   # ref
            set_allele.extend(re.split(',',fields_WGSvcf[4])) # alt

            if set_allele.index(major_allele) == 0: 
                ind_major_ref = 1; # indicator that ref. allele is major
            else:
                ind_major_ref = 0

            target_rv = {}
            for gt in gt_allele[0:2]:
                if gt == '.': continue
                elif set_allele[int(gt)] in rv_allele: # 
                    if len(target_rv) == 0:
                        target_rv[set_allele[int(gt)]] = 1
                    elif (len(target_rv) == 1) and (target_rv.keys()[0] == set_allele[int(gt)]):
                        target_rv[set_allele[int(gt)]] += 1
                    elif (len(target_rv) == 1) and (target_rv.keys()[0] != set_allele[int(gt)]):    
                        target_rv[set_allele[int(gt)]] = 1
            
            for ind_rv,count_var in target_rv.iteritems():
                fields_gen1k_var = []; out_site = []
                for gen1k_var in EURTabix.fetch(chrom,pos-1,pos):
                    fields_gen1k_var = gen1k_var.rstrip().split('\t')
                
                if len(fields_gen1k_var) == 0: # no variant in EUR population
                    out_site = [str(chrom),str(pos),str(fields_WGSvcf[3]),str(ind_rv),str(gencode_idx),str(gene_name), \
                            str(ind_major_ref),str(count_var),str(abs(int(tss_pos)-int(fields_WGSvcf[1])))]
                elif len(fields_gen1k_var) > 0:
                    for gen1k_allele in fields_gen1k_var[4:]:
                        temp_allele = gen1k_allele.split(':')
                        if (temp_allele[0] == ind_rv) and (temp_allele[1] < 0.01):
                            out_site = [str(chrom),str(pos),str(fields_WGSvcf[3]),str(ind_rv),str(gencode_idx),str(gene_name), \
                                    str(ind_major_ref),str(count_var),str(abs(int(tss_pos)-int(fields_WGSvcf[1])))]
                if len(out_site) > 0:
                    site_out.write("\t".join(out_site)+"\n")
        elif len(rv_allele) == 0:
            continue

#print("--- %s hours ---" % ((time.time() - start_time)/3600))

#f = open('${RAREVARDIR}/RIVER/data/score/indiv/done.txt', 'a+')
#f.write("%s\n" % (str(options.num)))
#f.close()
