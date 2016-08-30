#!/usr/bin/env python

## pick single tissue outliers by parsing flat file
## for each tissue-gene pair, this is the individual with the most extreme expression (and |Z| >= 2)
## for each tissue also make file where the most extreme individual is picked with no threshold requirement
## also for each tissue make file with matrix of which individuals were tested for which genes
## these last two sets of files will be used by the control picking script to make features

## outputs outliers as tab-delimited file with the following columns:
## GENE INDS TISSUE Z

infile = '../preprocessing/gtex_2015-01-12_normalized_expression.txt'
outfile = '../data/outliers_singlez_picked.txt'
noThreshPrefix = '../data/singlez/outliers_singlez_nothreshold_'
genefile = '../reference/gencode.v19.genes.v6p.patched_contigs_genetypes_autosomal.txt'
indivfile = '../data/outliers_medz_picked_counts_per_ind.txt'

zscores = open(infile, 'r')
outliers = open(outfile, 'w')
genetypes = open(genefile, 'r')
indivCounts = open(indivfile, 'r')
nothresh = dict() # for file handles (tuple); make these on the fly


# get list of pc and linc RNA genes
print "Getting list of genes to keep..."
genesToKeep = []
for geneline in genetypes:
    geneline = geneline.strip().split()
    if geneline[1] == 'protein_coding' or geneline[1] == 'lincRNA':
        genesToKeep.append(geneline[0])
genetypes.close()

# get list of individual with fewer than 50 multi-tissue outliers
print "Getting list of individuals to keep..."
indivsToKeep = []
for indivline in indivCounts:
    indivline = indivline.strip().split()
    if int(indivline[1]) < 50:
        indivsToKeep.append(indivline[0])
indivCounts.close()

header = zscores.readline().strip().split()
individuals = header[2:]

# retrict to individuals in keep list
keepInd = [ind in indivsToKeep for ind in individuals] # also need to filter future lines
individuals = [ind for keep, ind in zip(keepInd,individuals) if keep]

outliers.write('GENE\tINDS\tTISSUE\tZ\n')

## function to create file handles for a new tissue
## adds them to the nothresh dict
## note this function relies on pre-existing variables
## modifies the nothresh dict directly
def addHandles(tissuename):
    outlierFile = noThreshPrefix + tissuename + "_picked.txt"
    countFile = noThreshPrefix + tissuename + "_counts.txt"
    out = open(outlierFile, 'w')
    count = open(countFile, 'w')
    out.write('GENE\tINDS\tTISSUE\tZ\n')
    count.write('GENE\t' + '\t'.join(individuals) + "\n")
    nothresh[tissuename] = (out, count)

# each line is a tissue-gene pair. iterate through them
print "Identifying single-tissue outliers..."
for line in zscores:
    line = line.strip().split()
    # if not a protein-coding or lincRNA gene, skip
    if line[1] not in genesToKeep:
        continue
    tissue = line[0]
    gene = line[1]
    line = [val for keep, val in zip(keepInd, line[2:]) if keep]
    # get list of individuals in that tissue and then replace NA by 0
    countstring = gene + '\t' + '\t'.join(['0' if i == "NA" else '1' for i in line]) + '\n'
    values = [0 if i == "NA" else float(i) for i in line]
    # get most extreme individual
    index, value = max(enumerate(values), key = lambda x: abs(x[1]))
    outstring = '\t'.join([gene, individuals[index], tissue, str(value)]) + '\n'
    # write to combined outlier file if the individual passes the Z threshold
    if abs(value) >= 2:
        outliers.write(outstring)
    # write to files specific for this tissue
    if tissue not in nothresh:
        addHandles(tissue)
    out, count = nothresh[tissue]
    out.write(outstring)
    count.write(countstring)

zscores.close()
outliers.close()
for handle1, handle2 in nothresh.values():
    handle1.close()
    handle2.close()

print "Done!"
