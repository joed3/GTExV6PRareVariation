#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

# Load required packages and functions
require(ggplot2)
require(cowplot)

source('enrichment.functions.R')
options(stringsAsFactors=FALSE)

#------------------- FUNCTIONS
## Function that runs the generic read_file function specified in enrichment.functions.R
## and also creates "any" features
read_file_mod = function(dir, mafs, types, suffix, names, prefix) {
    count = read_file(dir, mafs, types, suffix, names, prefix)
    count$any_variant = (count$n_variants > 0) + 0
    return(count)
}

## Function to return the subset of the counts data frame that includes individuals
## with fewer than the specified number of outliers
filter_num_outliers = function(countdata, outlierdata, threshold) {
    exclude.inds = outlierdata$ind[outlierdata$noutliers >= threshold]
    exclude.genes = countdata$gene_id[countdata$ID %in% exclude.inds & countdata$Y == 1]
    counts.subset = countdata[!(countdata$ID %in% exclude.inds),]
    counts.subset = counts.subset[!(counts.subset$gene_id %in% exclude.genes), ]
    return(counts.subset)
}

## Function that runs proportion.ratios from enrichment.functions.R
## and also sets factor levels
counts2props = function(countdata, mafs, typenames) {
    prop.ratio = proportion.ratios(countdata)
    prop.ratio$MAF = factor(prop.ratio$MAF, levels = mafs)
    prop.ratio$TYPE = factor(prop.ratio$TYPE, levels = typenames)
    return(prop.ratio)
}

## Function to plot count enrichments by MAF facetted by variant type
## Include error bars
ratio.plot = function(plotdata, hasAlpha = FALSE, legend.title = '') {
    p = ggplot(plotdata, aes(x = MAF, y = ESTIM)) +
        theme_bw() + xlab('Minor allele frequency (%)') + ylab('Enrichment') +
        geom_abline(intercept = 1, slope = 0) +
        scale_colour_manual(values=colors.type) + guides(colour = FALSE) +
        facet_wrap( ~ TYPE, scales = "free") + 
        theme(strip.background = element_blank())
    if (hasAlpha) {
        p = p + geom_pointrange(aes(x = MAF, ymin = CI.LOW, ymax = CI.HIGH, colour = TYPE,
                                    alpha = GROUP, group = interaction(MAF, GROUP)),
                                position = position_dodge(width = 0.6)) +
            scale_alpha_discrete(range = c(1, 0.4), name = legend.title) +
            theme(legend.key = element_blank())
    } else {
        p = p + geom_pointrange(aes(x = MAF, ymin = CI.LOW, ymax = CI.HIGH, colour = TYPE))
    }
}



#------------------- MAIN

## Set variety of variables that will be useful
MAFs = c('0-1', '1-5', '5-10', '10-25')
TYPEs = c('SNPs', 'indels', 'HallLabSV')
type.names = c('SNV','Indel','SV')
colors.type = c('dodgerblue3','springgreen4','orangered3')
names(colors.type) = type.names

MAFs.top0 = paste0(MAFs, '.peer.top0')
MAFs.top5 = paste0(MAFs, '.peer.top5')
peer.types = c('Fully corrected', 'Top 5 PEER factors removed', 'No PEER factors removed')

removed.thresholds = c(500, 50, 30)
removed.names = c('None','outlier for >= 50 genes','outlier for >= 30 genes')

## Read in features and counts
OUT.DIR = paste0(dir, '/data/medz/')
# for fig 2b.
medz.count = read_file_mod(OUT.DIR, MAFs, TYPEs, '_counts_MAF', type.names, prefix = '10kb_')
# for supp. fig without PC
medz.count.withoutPC = read_file_mod(OUT.DIR, MAFs, TYPEs, '_counts_MAF', type.names, prefix = '10kb_noPC_')
# for supp. fig with different levels of peer correction
medz.count.top0 = read_file_mod(OUT.DIR, MAFs.top0, TYPEs, '_counts_MAF', type.names, prefix = '10kb_')
medz.count.top0$MAF = gsub('\\.peer\\.top0', '', medz.count.top0$MAF)
medz.count.top5 = read_file_mod(OUT.DIR, MAFs.top5, TYPEs, '_counts_MAF', type.names, prefix = '10kb_')
medz.count.top5$MAF = gsub('\\.peer\\.top5', '', medz.count.top5$MAF)

# for supp. fig with different number of individuals removed due to excess outliers
OUT.DIR.ALL = paste0(dir, '/data/medzall/')
medz.count.allinds = read_file_mod(OUT.DIR.ALL, MAFs, TYPEs, '_counts_MAF', type.names, prefix = '10kb_')
outlier.counts = read.table(paste0(dir, '/data/outliers_medz_picked_counts_per_ind.txt'), header = F,
                            col.names = c('ind', 'noutliers'))
medz.count.max30 = filter_num_outliers(medz.count.allinds, outlier.counts, 30)
medz.count.nomax = filter_num_outliers(medz.count.allinds, outlier.counts, 500)


## Get proportion ratios for each of the count groups above
medz.count.prop.ratio = counts2props(medz.count, MAFs, type.names)
medz.count.prop.ratio.withoutPC = counts2props(medz.count.withoutPC, MAFs, type.names)
medz.count.prop.ratio.top0 = counts2props(medz.count.top0, MAFs, type.names)
medz.count.prop.ratio.top5 = counts2props(medz.count.top5, MAFs, type.names)
medz.count.prop.ratio.max30 = counts2props(medz.count.max30, MAFs, type.names)
medz.count.prop.ratio.nomax = counts2props(medz.count.nomax, MAFs, type.names)


## MAIN FIGURE
count.ratio.plot = ratio.plot(medz.count.prop.ratio)

## SUPPLEMENTAL PLOT excluding individuals with more than 60,000 variants
varCounts = read.table(paste0(dir, '/data/variant_counts_per_individual_all.txt'), header = T, row.names = 1)
indIncl = rownames(varCounts)[rowSums(varCounts, na.rm = T) <= 60000]
medz.count.subset = medz.count[medz.count$ID %in% indIncl,]
medz.count.prop.ratio.subset = counts2props(medz.count.subset, MAFs, type.names)
medz.count.prop.ratio.subset$GROUP = "Excluding individuals\nwith > 60000 variants"
medz.count.prop.ratio$GROUP = "All individuals"
medz.count.prop.ratio.subset = rbind(medz.count.prop.ratio.subset, medz.count.prop.ratio)

count.ratio.plot.subset = ratio.plot(medz.count.prop.ratio.subset, hasAlpha = T)

## SUPPLEMENTAL PLOT with/without PC
medz.count.prop.ratio.withoutPC$GROUP = "Excluding\nexon variants"
medz.count.prop.ratio$GROUP = "All variants"
medz.count.prop.ratio.with.without = rbind(medz.count.prop.ratio, medz.count.prop.ratio.withoutPC)
medz.count.prop.ratio.with.without$GROUP = factor(medz.count.prop.ratio.with.without$GROUP)

count.ratio.plot.with.withoutPC = ratio.plot(medz.count.prop.ratio.with.without, hasAlpha = T)

## SUPPLEMENTAL PLOT different PEER correction levels
# adding peer column and combining peer ratios
medz.count.prop.ratio.top0$GROUP = factor(peer.types[3], levels = peer.types)
medz.count.prop.ratio.top5$GROUP = factor(peer.types[2], levels = peer.types)
medz.count.prop.ratio.topall = medz.count.prop.ratio
medz.count.prop.ratio.topall$GROUP = factor(peer.types[1], levels = peer.types)
medz.count.prop.ratio.peer = rbind(medz.count.prop.ratio.top0,
                              medz.count.prop.ratio.top5,
                              medz.count.prop.ratio.topall)

count.ratio.plot.peer = ratio.plot(medz.count.prop.ratio.peer, hasAlpha = T)

## SUPPLEMENTAL PLOT excluding different number of individuals
medz.count.prop.ratio.max30$GROUP = factor(removed.names[3], levels = removed.names)
medz.count.prop.ratio.nomax$GROUP = factor(removed.names[1], levels = removed.names)
medz.count.prop.ratio.max50 = medz.count.prop.ratio
medz.count.prop.ratio.max50$GROUP = factor(removed.names[2], levels = removed.names)
medz.count.prop.ratio.maxoutlier = rbind(medz.count.prop.ratio.max30,
                                         medz.count.prop.ratio.max50,
                                         medz.count.prop.ratio.nomax)

count.ratio.plot.outliers.thresh = ratio.plot(medz.count.prop.ratio.maxoutlier,
                                              hasAlpha = T, legend.title = 'Excluded individuals')

# Save workspace image
save.image(paste0(dir, '/data/figure2a.count.enrichments.RData'))

fsize = 8
count.ratio.plot.peer = count.ratio.plot.peer +
    theme(axis.text = element_text(size = fsize),
          axis.title = element_text(size = fsize),
          legend.text = element_text(size = fsize),
          strip.text = element_text(size = fsize + 1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())

labeled.peer = ggdraw() +
    draw_plot(count.ratio.plot.peer, 0,0,1,1) +
    draw_plot_label('e', 0, 0.97, size = fsize + 3)

pdf(paste0(dir, '/paper_figures/suppfig.peer.enrichments.panel.pdf'), width = 9, height = 2.6)
labeled.peer
dev.off()
