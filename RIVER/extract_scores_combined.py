#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author - Yungil Kim
:E-mail - ipw012@gmail.com

:Description:

From each of individual.txt files, This script extract all of features including chromHMM, Phylop, DANN scores and combine all of features simultaneously

cat ${RAREVARDIR}/RIVER/data/score/indiv/${ID}.txt |
${RAREVARDIR}/RIVER/code/extract_scores_combined.py -n $count_ind --id $ID
--af_in ${RAREVARDIR}/data/wgs/GTEx_af.vcf.gz --wgs_in
filtered_and_compressed_GTEx_WGS_vcf_file --anno_in
${RAREVARDIR}/data/wgs/GTEx_vep.vcf.gz --cadd_in
${RAREVARDIR}/RIVER/data/CADD_whole_genome_SNVs_inclAnno.tsv.gz --dann_in
${RAREVARDIR}/RIVER/data/DANN_whole_genome_SNVs.tsv.bgz --chromHMM_in
${RAREVARDIR}/RIVER/data/wgEncodeBroadHmmGm12878HMM.sorted.hg19.bed.txt.gz
--phylop_in ${RAREVARDIR}/RIVER/data/phyloP100way.txt.gz --score_out
${RAREVARDIR}/RIVER/data/score/indiv/${ID}.${count_ind}.score.nuc.txt

"""

import sys, os, re
import pysam
import numpy as np
from optparse import OptionParser
import time

# start_time = time.time()

parser = OptionParser()

parser.add_option("-n", type="int", dest="num")
parser.add_option("--id",action="store", type="string", dest="ind_id")

parser.add_option("--af_in", dest="af_in")          # af vcf GTEX
parser.add_option("--wgs_in", dest="wgs_in")        # wgs vcf
parser.add_option("--anno_in", dest="anno_in")    # variant annotation

parser.add_option("--cadd_in", dest="cadd_in")    # CADD.tsv
parser.add_option("--dann_in", dest="dann_in")          # DANN core
parser.add_option("--chromHMM_in", dest="chromHMM_in")  # chromHMM
parser.add_option("--phylop_in", dest="phylop_in")      # phylop

# output
parser.add_option("--score_out", dest="score_out")  # Deleteriousness

(options, args) = parser.parse_args()

if os.path.exists(options.wgs_in) and os.path.exists(options.wgs_in+".tbi"):
    wgsTabix = pysam.Tabixfile(options.wgs_in,'r')                  # GTEx VCF file
if os.path.exists(options.af_in) and os.path.exists(options.af_in+".tbi"):
    afTabix = pysam.Tabixfile(options.af_in,'r')                    # GTEx AF file
if os.path.exists(options.anno_in) and os.path.exists(options.anno_in+".tbi"):
    annoTabix = pysam.Tabixfile(options.anno_in,'r')                    # variant annotation
if os.path.exists(options.cadd_in) and os.path.exists(options.cadd_in+".tbi"):
    caddTabix = pysam.Tabixfile(options.cadd_in,'r')                    # CADD.tsv
if os.path.exists(options.dann_in) and os.path.exists(options.dann_in+".tbi"):
    dannTabix = pysam.Tabixfile(options.dann_in,'r')                    # DANN score file
if os.path.exists(options.chromHMM_in) and os.path.exists(options.chromHMM_in+".tbi"):
    chromHMMTabix = pysam.Tabixfile(options.chromHMM_in,'r')            # chromHMM
if os.path.exists(options.phylop_in) and os.path.exists(options.phylop_in+".tbi"):
    phylopTabix = pysam.Tabixfile(options.phylop_in,'r')                # phylop

WGS_header = []
with open(options.wgs_in,"r") as f:
    for line in f:
        eachcol = line.rstrip().split('\t')
        if eachcol[0].startswith("#CHROM") == 1:
            for colnames in eachcol:
                WGS_header.append(colnames)
        break

dic_wgs = {}
for col_vcfs in WGS_header:
    if len(col_vcfs) > 9:
        dic_wgs[col_vcfs[0:9]] = WGS_header.index(col_vcfs)
    else:
        dic_wgs[col_vcfs] = WGS_header.index(col_vcfs)

dic_chromHMM = {"1_Active_Promoter": 1, "2_Weak_Promoter": 2, "3_Poised_Promoter": 3, "4_Strong_Enhancer": 4, \
                "5_Strong_Enhancer": 4, "6_Weak_Enhancer": 5, "7_Weak_Enhancer": 5, "8_Insulator": 6, \
                "9_Txn_Transition": 7, "10_Txn_Elongation": 7, "11_Weak_Txn": 8, "12_Repressed": 9, \
                "13_Heterochrom/lo": 10, "14_Repetitive/CNV": 10, "15_Repetitive/CNV": 10}

dic_segway = {"C0": 1, "C1": 1, "D": 2, "E/GM": 3, "F0": 4, "F1": 4, "GE0": 5, "GE1": 5, "GE2": 5, "GM0": 6, "GM1": 6, \
              "GS": 7, "H3K9me1": 8, "L0": 9, "L1": 9, "R0": 10, "R1": 10, "R2": 10, "R3": 10, "R4": 10, "R5": 10, \
              "TF0": 11, "TF1": 11, "TF2": 11, "TSS": 12}

# Complete header from CADD data
cadd_header = []
with open(options.cadd_in,"r") as f:
    for line in f:
        eachcol = line.rstrip().split('\t')
        if eachcol[0].startswith("#Chrom") == 1:
            for colnames in eachcol:
                cadd_header.append(colnames)
        break

score_header = ["Ensembl_id",'anno',"Chrom","Pos","nvar",'GC','CpG','priPhCons','mamPhCons','verPhCons','priPhyloP', \
                'mamPhyloP','verPhyloP','GerpN','GerpS','dnaHelT','dnaMGW','dnaProT','dnaRoll','fitCons', \
                'cHmmTssA','cHmmTssAFlnk','cHmmTxFlnk','cHmmTx','cHmmTxWk','cHmmEnhG','cHmmEnh','cHmmZnfRpts','cHmmHet', \
                'cHmmTssBiv','cHmmBivFlnk','cHmmEnhBiv','cHmmReprPC','cHmmReprPCWk','cHmmQuies','EncH3K27Ac','EncH3K4Me1', \
                'EncH3K4Me3','EncNucleo','EncOCCombPVal','EncOCDNasePVal','EncOCFairePVal','EncOCpolIIPVal','EncOCctcfPVal', \
                'EncOCmycPVal','EncOCDNaseSig','EncOCFaireSig','EncOCpolIISig','EncOCctcfSig','EncOCmycSig','TFBS','TFBSPeaks', \
                'TFBSPeaksMax','PHRED','DistTSS','Segway',"chromHMM","phylop","DANN"]

# idx_labels = []
# for features in score_header[5:-3]:
#     if features == "DistTSS": continue # Distance to TSS is going to be added at the final step
#     else: idx_labels.extend([CADDtsv_header.index(features)])

score_out = open(options.score_out,'w');  score_out.write("\t".join(score_header)+"\n")

# read lines from text files for reading targeted regions
# chr | pos | ref | rv_allele | gencode_idx | ensembl_id | gene name | ind_major_ref | count_rv | distTSS
# ensembl_id | chr | pos
for target_region in sys.stdin:
    fields_site = target_region.rstrip().split('\t')

    dic_header = {}

    dic_header["Chrom"] = fields_site[1]
    dic_header["Pos"] = fields_site[2]
    # ref = fields_site[2]
    # rv_allele = fields_site[3]
    # dic_header["GencodeIdx"] = fields_site[4]
    dic_header["Ensembl_id"] = fields_site[0]
    # dic_header["GeneName"] = fields_site[6]
    # ind_ref_major = fields_site[7] # indicator of reference allele == major
    # dic_header["nvar"] = fields_site[8]
    # dic_header["DistTSS"] = fields_site[9]; count_feature = 1

    # background should be also n_var*background
    # out_score = [str(gene_name), str(idx_gencode), str(chrom), str(pos), str(n_var)]

    dic_header["nvar"] = "NA"
    for poss_var in wgsTabix.fetch(dic_header["Chrom"],int(dic_header["Pos"])-1,int(dic_header["Pos"])):
        fields_WGSvcf = poss_var.rstrip().split('\t')

        ref = fields_WGSvcf[3]

        set_allele = [fields_WGSvcf[3]]   # ref
        set_allele.extend(re.split(',',fields_WGSvcf[4])) # alt

        for target_site in afTabix.fetch(dic_header["Chrom"],int(dic_header["Pos"])-1,int(dic_header["Pos"])):
            fields_af = target_site.rstrip().split('\t')

            minor_allele = []
            for gtex_allele in fields_af[4:]:
                allele_af = gtex_allele.split(':')
                if (len(allele_af[0]) == 1) and (allele_af[0] == ref) and (float(allele_af[1]) > 0.01): # SNV & reference major
                    ind_ref_major = 1
                elif (len(allele_af[0]) == 1) and (allele_af[0] == ref) and (float(allele_af[1]) <= 0.01) and (float(allele_af[1]) > 0): # SNV & reference minor
                    ind_ref_major = 0
                if (float(allele_af[1]) <= 0.01) and (float(allele_af[1]) > 0):
                    minor_allele.extend([allele_af[0]])

        gt_ind = fields_WGSvcf[int(dic_wgs[options.ind_id])] # individual-specific allele
        gt_allele = re.split('[/:]',gt_ind)
        # print '%s' %(gt_ind)
#       print '%s' %(minor_allele)
        count_var = 0
        for ma in minor_allele:
            # print '%s' %(ma)
            for gt in gt_allele[0:2]:
                if gt == '.': continue
                elif int(gt) == set_allele.index(ma): rv_allele = set_allele[int(gt)]; count_var += 1
        # print '%d' %(count_var)
        dic_header["nvar"] = count_var

    list_tss = open('${RAREVARDIR}/reference/gencode.v19.genes.v6p.patched_contigs_TSS.bed','r')
    for tss_pos in list_tss.readlines():
        tss_info = tss_pos.rstrip().split('\t')
        if tss_info[3] == dic_header["Ensembl_id"]:
            tss = tss_info[2]

    dic_header["DistTSS"] = abs(int(tss)-int(dic_header["Pos"])) ; count_feature = 1

    # possible site from CADD data
    if int(ind_ref_major) == 1: # ref is major allele
        for anno_site in annoTabix.fetch(dic_header["Chrom"],int(dic_header["Pos"])-1,int(dic_header["Pos"])):
            fields_anno = anno_site.rstrip().split('\t')
            temp_field1 = re.sub("CSQ=",",",fields_anno[-1])
            temp_field2 = re.split('[,]',temp_field1) # per gene

            list_anno = [];
            for fields_bygene in temp_field2[1:]:
                temp_field3 = re.split('[|]',fields_bygene)
                if temp_field3[4] == re.split('[.]',dic_header["Ensembl_id"])[0]: # same gene
                    list_anno.extend([temp_field3[1]]);
                else: continue

            count_feature += 1
            if len(list_anno) > 0:
                dic_header["anno"] = "&".join(list_anno);
            else:
                dic_header["anno"] = "undefined"

        for poss_site in caddTabix.fetch(dic_header["Chrom"],int(dic_header["Pos"])-1,int(dic_header["Pos"])):
            #if poss_var.startswith('#'): continue
            fields_cadd = poss_site.rstrip().split('\t')

            if (ref == fields_cadd[cadd_header.index("Ref")]) and (rv_allele == fields_cadd[cadd_header.index("Alt")]):   # match ref and alt allele
                for feature in score_header[5:-4]:
                    if feature == "DistTSS": continue
                    elif feature == "PHRED":
                        dic_header[feature] = float(fields_cadd[cadd_header.index(feature)]); count_feature += 1;
                    elif fields_cadd[cadd_header.index(feature)] == "NA":
                        dic_header[feature] = "NaN"
                    elif feature == "Segway":
                        dic_header[feature] = dic_segway[str(fields_cadd[cadd_header.index(feature)])]; count_feature += 1
                    else:
                        dic_header[feature] = float(fields_cadd[cadd_header.index(feature)]); count_feature += 1
                break
            else: continue
        # DANN
        dic_header["DANN"] = "NaN"
        for dann_site in dannTabix.fetch(dic_header["Chrom"],int(dic_header["Pos"])-1,int(dic_header["Pos"])):
            fields_dann = dann_site.rstrip().split('\t')
            if (fields_dann[2] == ref) and (fields_dann[3] == rv_allele):
                dic_header["DANN"] = float(fields_dann[4]); count_feature += 1

    elif int(ind_ref_major) == 0: # ref is minor allele
        dic_header["anno"] = "undefined"

        for poss_site in caddTabix.fetch(dic_header["Chrom"],int(dic_header["Pos"])-1,int(dic_header["Pos"])):
            #if poss_var.startswith('#'): continue
            fields_cadd = poss_site.rstrip().split('\t')
            if ref == fields_cadd[cadd_header.index("Ref")]:
                for feature in score_header[5:-4]:
                    if feature == "DistTSS": continue
                    elif feature == "PHRED":
                        dic_header[feature] = "NaN"
                    elif fields_cadd[cadd_header.index(feature)] == "NA":
                        dic_header[feature] = "NaN"
                    elif feature == "Segway":
                        dic_header[feature] = dic_segway[str(fields_cadd[cadd_header.index(feature)])]; count_feature += 1
                    else:
                        dic_header[feature] = float(fields_cadd[cadd_header.index(feature)]); count_feature += 1
                break
        # DANN
        dic_header["DANN"] = "NaN"

    # chromHMM
    dic_header["chromHMM"] = "NaN"
    for poss_site in chromHMMTabix.fetch(dic_header["Chrom"],int(dic_header["Pos"])-1,int(dic_header["Pos"])):
        fields_chromHMM = poss_site.rstrip().split('\t')
        dic_header["chromHMM"] = dic_chromHMM[str(fields_chromHMM[11])]; count_feature += 1

    # phylop
    dic_header["phylop"] = "NaN"
    for poss_site in phylopTabix.fetch(dic_header["Chrom"],int(dic_header["Pos"])-1,int(dic_header["Pos"])):
        fields_phylop = poss_site.rstrip().split('\t')
        dic_header["phylop"] = float(fields_phylop[2]); count_feature += 1

    if count_feature > 1:
        out_score = []
        for feature in score_header:
            out_score.extend([str(dic_header[feature])])
        score_out.write("\t".join(out_score)+"\n")
#print("--- %s hours ---" % ((time.time() - start_time)/3600))

#f = open('${RAREVARDIR}/RIVER/data/score/done.combined.txt', 'a+')
#f.write("%s\n" % (str(options.num)))
#f.close()
