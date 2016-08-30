#!/use/bin/evn Rscript

library(data.table)
library(ggplot2)

# Look into into overlap of single-tissue outliers between tissue
# and into the overlap of single-tissue and multi-tissue outliers
# Also look at number of single-tissue outliers per tissue, per individual

medz = fread('../data/outliers_medz_picked.txt', header = T, stringsAsFactors = F)
singlez = fread('../data/outliers_singlez_picked.txt', header = T, stringsAsFactors = F)
setkey(singlez, GENE)

# number of multi-tissue outliers per individual
cat("Median number of genes for which an individual is a multi-tissue outliers",
    median(table(medz$INDS)), "\n")

# overlap of single-tissue and multi-tissue outliers
setkey(medz, GENE, INDS)
setkey(singlez, GENE, INDS)

single.med = merge(singlez, medz, all.x = T)
overlap = single.med[, .(total = nrow(.SD), multi = nrow(.SD[!is.na(DFS),])), by = .(TISSUE)]
overlap$prop = overlap$multi / overlap$total

cat("minimum overlap between single-tissue and multi-tissue outliers: ", min(overlap$prop), "\n")
cat("maximum overlap between single-tissue and multi-tissue outliers: ", max(overlap$prop), "\n")

# read in the individual by tissue associations
tissue.ind = fread('../gtex_2015-01-12_tissue_by_ind.txt', header = T, stringsAsFactors = F)
# not removing the individuals with too many outliers because they were included in the Z-score calculations...
samples.tissue = tissue.ind[, .(nsamples = nrow(.SD)), by = Tissue]

# plot relationship between tissue sample size and overlap with median Z
setnames(samples.tissue, c('Tissue'), c('TISSUE'))
setkey(samples.tissue, TISSUE)
setkey(overlap, TISSUE)

samples.tissue.overlap = merge(samples.tissue, overlap)

correlation = cor.test(samples.tissue.overlap$nsamples, samples.tissue.overlap$prop)
print(correlation)

# some stats on the single-tissue outliers
cat('number of outliers per tissue:')
print(overlap[order(total), .(TISSUE, total)])

cat('median number of outliers per tissue:', median(overlap$total), '\n')

singlez.per.tissue.ind = singlez[, .(noutliers = nrow(.SD)), by = .(TISSUE, INDS)]
singlez.per.tissue = singlez.per.tissue.ind[, .(total = sum(.SD$noutliers), median = as.double(median(.SD$noutliers))), by = INDS]
cat('median number of outliers per individual:', median(singlez.per.tissue$total), '\n')
cat('median number of outliers per tissue per individual:', median(singlez.per.tissue$median), '\n')

# get number of outliers as a percentage of tested genes
overlap$single.tested = integer(44)
for (i in 1:nrow(overlap)) {
    fname = paste0('../singlez/outliers_singlez_nothreshold_',
        overlap$TISSUE[i], '_picked.txt')
    overlap[i, 'single.tested'] = nrow(fread(fname))
}

cat('percentage of genes that have outliers per tissue:')
print(overlap[order(total/single.tested), total/single.tested])

save.image('data/suppfig.single.multi.overlap.RData')

