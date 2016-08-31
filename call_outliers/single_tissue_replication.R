#!/usr/bin/env Rscript

library(data.table)

### Master directory
dir = Sys.getenv('RAREVARDIR')

set.seed(1)
# test whether replication of single-tissue outliers is driven by shared individuals among closely shared tissues
# here, we define replication as in figure 1b:
# the outlier in tissue 1 (most extreme individual, requiring |Z| >= 2) has |Z| >= 2 in tissue 2 and in same direction

# based on a small simulation, it looks like it is ok to stardardize subsets of already standardized values
# leads to the same result as subsetting than standardizing


# FUNCTIONS ------------------------------
# takes a data frame standardizes the rows (to have mean 0 and variance 1)
# returns the transformed data frame
standardize = function(data) {
    transformed = t(apply(data, 1, function(x) (x - mean(x)) / sd(x)))
    return(as.data.frame(transformed))
}

# takes in data frames of expr values
# returns replication fo t1 in t2 (discover in t1, replicate in t2)
# requires that t1 and t2 have the same individuals (columns) in the same order
replication.shared = function(t1, t2, thresh = 2) {
    stopifnot(ncol(t1) == ncol(t2))
    stopifnot(sum(colnames(t1) == colnames(t2)) == ncol(t1))
    t1 = standardize(t1)
    t2 = standardize(t2)
    t1.indx = apply(abs(t1), 1, which.max)
    # check replication in tissue 2
    t2.repli.z = mapply(function(i,j) t2[i,j] * sign(t1[i,j]), c(1:(length(t1.indx))), t1.indx)
    # only considering replication for cases where the |Z-score| in tissue 1 was >= 2
    keep = mapply(function(i,j) abs(t1[i,j]) >= 2, c(1:(length(t1.indx))), t1.indx)
    replication = sum(t2.repli.z[keep] >= thresh) / sum(keep)
    return(replication)
}

# takes unstandardized data frames of expr values from t1 and t2 (shared individuals between them)
# as well as the unstandardized expr values of t2 (a distinct set of individuals)
# returns the replication of the outliers identified in t1 in the distinct t2 individuals
replication.disjoint = function(t1.shared, t2.shared, t2, thresh = 2) {
    t1.indx = apply(abs(standardize(t1.shared)), 1, which.max) # get outliers in discovery tissue
    t2.values = mapply(function(i,j) t2.shared[i,j], c(1:(length(t1.indx))), t1.indx) # get outlier values from shared in t2
    t2$t1.outlier = t2.values # append these to t2
    t2 = standardize(t2) # standardize t2 and check whether outliers have |Z-score| >= 2
    # only consider replication for cases with |Z-score| >= 2 in tissue 1
    keep = mapply(function(i,j) abs(t1.shared[i,j]) >= 2, c(1:(length(t1.indx))), t1.indx)
    replication = sum(t2$t1.outlier[keep] >= thresh) / sum(keep)
    return(replication)
}

# get the maximum group size such that we can have n-1 in each of t1 only and t2 only
# and n in the shared group
# return the maximum n, where n is one of 125, 100, 70
# if the max is less than 70, return NA
get.group.size = function(nshared, nt1, nt2, maxsize) {
    # limited by shared
    selected = NA
    if (nshared <= (min(nt1, nt2) + 1)) {
        if (nshared >= 125) {
            selected = 125
        } else if (nshared >= 100) {
            selected = 100
        } else if (nshared >= 70) {
            selected = 70
        }
    # shared > smallest of either tissue
    } else {
        # loop through different sizes from largest to smallest
        for (n in c(125, 100, 70)) {
            if (nshared < n) next
            excess = nshared - n
            missing = max(0, (n - 1 - nt1)) + max(0, (n - 1 - nt2))
            if (excess >= missing) {
                selected = n
                break
            }
        }
    }
    if (!is.na(selected) & selected > maxsize) {
        selected = maxsize
    }
    return(selected)
}

# split individuals into those that will be in the shared group
# and those into each of the "tissue-only" groups
# note that some of the shared individuals will be split between the tissue-only groups
get.groups = function(tis1.ind, tis2.ind, max.group.size) {
    shared = intersect(tis1.ind, tis2.ind)
    tis1only = setdiff(tis1.ind, tis2.ind)
    tis2only = setdiff(tis2.ind, tis1.ind)
    n = get.group.size(length(shared), length(tis1only), length(tis2only), maxsize = max.group.size)
    if (is.na(n)) {
        return(NULL)
    }
    # redistribute shared individuals to maximize the group sizes
    # sample down to n-1 in each group (n for shared)
    # shuffle shared then split it up
    shared = sample(shared)
    tis1.missing.n = n - length(tis1only) - 1
    tis2.missing.n = n - length(tis2only) - 1
    if (tis1.missing.n > 0) {
        tis1only = c(tis1only, shared[1:tis1.missing.n])
        shared = shared[-c(1:tis1.missing.n)]
    } else {
        tis1only = tis1only[1:(n-1)]
    }
    if (tis2.missing.n > 0) {
        tis2only = c(tis2only, shared[1:tis2.missing.n])
        shared = shared[-c(1:tis2.missing.n)]
    } else {
        tis2only = tis2only[1:(n-1)]
    }
    shared = shared[1:n]
    # some sanity checks
    stopifnot(length(tis1only) == (n-1))
    stopifnot(length(tis2only) == (n-1))
    stopifnot(length(shared) == n)
    stopifnot(length(intersect(tis1only, tis2only)) == 0)
    stopifnot(length(intersect(tis1only, shared)) == 0)
    stopifnot(length(intersect(tis2only, shared)) == 0)
    return(list(n = n, tis1only = tis1only, tis2only = tis2only, shared = shared))
}

# first identify outliers in the shared set of individuals in one tissue
# test replication in shared set in the other tissue
# and in distinct set of individuals in both tissues (unless shared.only = TRUE)
#
# returns a vectors with the number of samples in the comparison
# and the various replication rates
# (in order)
# discovery replication
# ---------------------
#    t1         t2 (same individuals)
#    t2         t1 (same individuals)
#    t1         t2
#    t1         t1 
#    t2         t1
#    t2         t2
# only returns the first 2 replication rates if shared.only = TRUE
run.replication = function(tissue.names, max.group.size = 125, shared.only = FALSE) {
    tname1 = tissue.names[1]
    tname2 = tissue.names[2]
    tis1 = tissue.expr[[tname1]]
    tis2 = tissue.expr[[tname2]]
    # remove genes that aren't in both tissues and order by gene name
    genes.both = c(tis1$Gene, tis2$Gene)[duplicated(c(tis1$Gene, tis2$Gene))]
    setkey(tis1, Gene) # setting the key takes care of the sorting
    setkey(tis2, Gene)
    tis1 = tis1[genes.both, ]
    tis2 = tis2[genes.both, ]
    stopifnot(nrow(tis1) == nrow(tis2))
    stopifnot(sum(tis1$Gene == tis2$Gene) == nrow(tis1))
    # split into shared and disjoint set of individuals
    tis1.ind = colnames(tis1)[-c(1,2)]
    tis2.ind = colnames(tis2)[-c(1,2)]
    if (shared.only) {
        inds = sample(intersect(tis1.ind, tis2.ind))
        if (length(inds) < max.group.size) {
            return(rep(NA,5))
        }
        inds = inds[1:max.group.size]
        tis1.shared.expr = as.data.frame(tis1[, inds, with = F])
        tis2.shared.expr = as.data.frame(tis2[, inds, with = F])
        # discover outliers in shared and test replication in shared in the other tissue
        tis1.in.tis2.shared = replication.shared(tis1.shared.expr, tis2.shared.expr)
        tis2.in.tis1.shared = replication.shared(tis2.shared.expr, tis1.shared.expr)
        return(c(tname1, tname2, max.group.size, tis1.in.tis2.shared, tis2.in.tis1.shared))
    } else {
        groups = get.groups(tis1.ind, tis2.ind, max.group.size)
        if (is.null(groups)) {
            return(rep(NA, 9))
        }
        # split expression data into the shared and not shared for each tissue
        tis1.shared.expr = as.data.frame(tis1[, c(groups$shared), with = F])
        tis2.shared.expr = as.data.frame(tis2[, c(groups$shared), with = F])
        tis1.only.expr = as.data.frame(tis1[, c(groups$tis1only), with = F])
        tis2.only.expr = as.data.frame(tis2[, c(groups$tis2only), with = F])
        # discover outliers in shared and test replication in shared in the other tissue
        tis1.in.tis2.shared = replication.shared(tis1.shared.expr, tis2.shared.expr)
        tis2.in.tis1.shared = replication.shared(tis2.shared.expr, tis1.shared.expr)
        # as well as in each of the non-shared groups
        tis1.in.tis2.only = replication.disjoint(tis1.shared.expr, tis2.shared.expr, tis2.only.expr)
        tis1.in.tis1.only = replication.disjoint(tis1.shared.expr, tis1.shared.expr, tis1.only.expr)
        tis2.in.tis1.only = replication.disjoint(tis2.shared.expr, tis1.shared.expr, tis1.only.expr)
        tis2.in.tis2.only = replication.disjoint(tis2.shared.expr, tis2.shared.expr, tis2.only.expr)
    
        return(c(tname1, tname2, groups$n, tis1.in.tis2.shared, tis2.in.tis1.shared,
                 tis1.in.tis2.only, tis1.in.tis1.only, tis2.in.tis1.only, tis2.in.tis2.only))
    }
}

# function to format the replication data frames nicely
# this function will probably need to be modified if start producing different types of df
# currently can handle the complete shared/not shared and shared only ones
clean.replication.df = function(df){
    column.names = c('tissue1', 'tissue2', 'n', 'tis1.in.tis2.shared', 'tis2.in.tis1.shared',
                     'tis1.in.tis2.only', 'tis1.in.tis1.only', 'tis2.in.tis1.only', 'tis2.in.tis2.only')
    column.names = column.names[1:ncol(df)] # for shared only dfs
    colnames(df) = column.names
    df = df[!is.na(df$n), ]
    for (i in c(3:ncol(df))) {
        df[,i] = as.numeric(df[, i])
    }
    return(df)
}

# ----------------------------------------

## read in flat file and subset to protein-coding genes
genetypes = read.table(paste0(dir, '/reference/gencode.v19.genes.v6p.patched_contigs_genetypes_autosomal.txt'),
    col.names = c('gene','type'))
keepGenes = genetypes$gene[genetypes$type %in% c('protein_coding','lincRNA')]
# read in expression data
expr = fread(paste0(dir, '/preprocessing/gtex_2015-01-12_normalized_expression.txt'))
expr = expr[Gene %in% keepGenes,]
# read in set of individuals that need to be exluded because of too many multi-tissue outliers
outlier.counts = read.table(paste0(dir, '/data/outliers_medz_picked_counts_per_ind.txt'), header = FALSE, stringsAsFactors = F)
keepInds = outlier.counts[outlier.counts[,2] < 50, 1]
expr = expr[, colnames(expr) %in% c('Tissue','Gene',keepInds), with = F]

## break flat file into a file for every tissue
tissues = unique(expr$Tissue)
tissue.expr = lapply(tissues, function(tis) expr[Tissue == tis, ])
names(tissue.expr) = tissues
# remove individuals from each tissue that have no data (based on first row)
# and count the number of individuals per tissue
tissue.counts = data.frame(tissue = character(length(tissues)), nsample = numeric(length(tissues)),
                           stringsAsFactors = FALSE)
i = 1
for (t in tissues) {
    dt = tissue.expr[[t]]
    tissue.expr[[t]] = dt[, !is.na(dt[1,]), with = FALSE]
    tissue.counts[i,1] = t
    tissue.counts[i,2] = ncol(tissue.expr[[t]]) - 2
    i = i + 1
}

# make list of pairs of tissues
pairs = t(combn(tissue.counts$tissue, 2))

# replication in 70 shared individuals
replications.shared70 = as.data.frame(t(apply(pairs, 1, run.replication, max.group.size = 70, shared.only = T)), stringsAsFactors = F)
# clean up data frames (remove NAs and make columns numeric)
replications.shared70 = clean.replication.df(replications.shared70)

# replication in shared and distinct individuals ** check/modify the get.group.size function if you change the max.group size
replications70 = as.data.frame(t(apply(pairs, 1, run.replication, max.group.size = 70)), stringsAsFactors = F)
replications70 = clean.replication.df(replications70)

save(replications.shared70, replications70, file = paste0(dir, '/data/single_tissue_replication.RData'))
