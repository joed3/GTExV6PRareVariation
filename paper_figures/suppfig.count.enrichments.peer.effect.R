#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

# Load required packages and functions
require(ggplot2)
require(gtable)
options(stringsAsFactors=FALSE)
#--------------------- FUNCTIONS

# for each variant type and each minor allele frequency bin,
# get the ratio of the portions of outliers/non outliers with rare variants
# and the associated confidence intervals (exp of Wald interval on the log of the proportion ratio)
proportion.ratios = function(counts) {
    mafs = unique(counts$MAF)
    types = unique(counts$TYPE)
    # make emty data frame and fill it from rbinds (not very many, so it's fine)
    results = data.frame(ESTIM=numeric(), CI.LOW=numeric(), CI.HIGH=numeric(), MAF=character(), TYPE=character(),
        stringsAsFactors = F)
    for (m in mafs) {
        for (t in types) {
            results = rbind(results, proportion.ratios.helper(counts, m, t))
        }
    }
    return(results)
}

# actually does the work of the main function for the given maf bin and variant type
proportion.ratios.helper = function(counts, maf, type) {
    counts.subset = counts[counts$MAF == maf & counts$TYPE ==type, c('Y','any_variant')]
    summary.counts = as.data.frame(table(counts.subset))
    stopifnot(nrow(summary.counts)== 4 & min(summary.counts$Freq) > 0) # check that no zero values, and right number of values
    # get required values to get statistic and CI
    out.var = summary.counts$Freq[summary.counts$Y == 1 & summary.counts$any_variant == 1]
    nonout.var = summary.counts$Freq[summary.counts$Y == 0 & summary.counts$any_variant == 1]
    out.total = sum(summary.counts$Freq[summary.counts$Y == 1])
    nonout.total = sum(summary.counts$Freq[summary.counts$Y == 0])
    # actually calculate
    estimate = (out.var/out.total)/(nonout.var/nonout.total)
    # get bounds of confidence interval on the log of the proportion then exponentiate
    log.se = sqrt(1/out.var - 1/out.total + 1/nonout.var - 1/nonout.total)
    max.ci = estimate * exp(1.96*log.se)
    min.ci = estimate * exp(-1.96*log.se)
    # put all together in a list that can become the row of a dataframe
    dfrow = list(ESTIM=estimate, CI.LOW=min.ci, CI.HIGH=max.ci, MAF=maf, TYPE=type)
    return (dfrow)
}

#------------------- MAIN

# Read in count data for fully corrected outliers/controls and those with 0/5 PEER factors removed
MAFs = c('0-1', '1-5', '5-10', '10-25')
MAFs.top0 = paste(MAFs, '.peer.top0', sep = '')
MAFs.top5 = paste(MAFs, '.peer.top5', sep = '')
TYPEs = c('SNPs', 'indels', 'HallLabSV')
type.names = c('SNV','Indel','SV')
peer.types = c('Fully corrected', 'Top 5 PEER factors removed', 'No PEER factors removed')
peer.alphas = c(1, .7, .4)
names(peer.alphas) = peer.types

OUT.DIR = paste0(dir, '/data/medz/')
medz.count = read_file(OUT.DIR, MAFs, TYPEs, '_counts_MAF', type.names, prefix = '10kb_')
medz.count$PEER = factor(peer.types[1], levels = peer.types)

medz.count.top0 = read_file(OUT.DIR, MAFs.top0, TYPEs, '_counts_MAF', type.names, prefix = '10kb_')
medz.count.top0$PEER = factor(peer.types[3], levels = peer.types)
medz.count.top0$MAF = gsub('\\.peer\\.top0', '', medz.count.top0$MAF)

medz.count.top5 = read_file(OUT.DIR, MAFs.top5, TYPEs, '_counts_MAF', type.names, prefix = '10kb_')
medz.count.top5$PEER = factor(peer.types[2], levels = peer.types)
medz.count.top5$MAF = gsub('\\.peer\\.top5', '', medz.count.top5$MAF)

# Turn counts into binary presence/absence 
medz.count$any_variant = (medz.count$n_variants > 0) + 0
medz.count.top0$any_variant = (medz.count.top0$n_variants > 0) + 0
medz.count.top5$any_variant = (medz.count.top5$n_variants > 0) + 0

# Make plot of count enrichments by MAF faceted by variant type
# Include error bars
colors.type = c('dodgerblue3','springgreen4','orangered3')
names(colors.type) = type.names

# Make plot of count enrichment for SNPs only using ratio of proporions instead of log odds
medz.count.prop.ratio = proportion.ratios(medz.count)
medz.count.prop.ratio$MAF = factor(medz.count.prop.ratio$MAF, levels = MAFs)
medz.count.prop.ratio$TYPE = factor(medz.count.prop.ratio$TYPE, levels = type.names)
medz.count.prop.ratio$PEER = factor(peer.types[1], levels = peer.types)

medz.count.prop.ratio.top0 = proportion.ratios(medz.count.top0)
medz.count.prop.ratio.top0$MAF = factor(medz.count.prop.ratio.top0$MAF, levels = MAFs)
medz.count.prop.ratio.top0$TYPE = factor(medz.count.prop.ratio.top0$TYPE, levels = type.names)
medz.count.prop.ratio.top0$PEER = factor(peer.types[3], levels = peer.types)

medz.count.prop.ratio.top5 = proportion.ratios(medz.count.top5)
medz.count.prop.ratio.top5$MAF = factor(medz.count.prop.ratio.top5$MAF, levels = MAFs)
medz.count.prop.ratio.top5$TYPE = factor(medz.count.prop.ratio.top5$TYPE, levels = type.names)
medz.count.prop.ratio.top5$PEER = factor(peer.types[2], levels = peer.types)

medz.count.prop.ratio = rbind(medz.count.prop.ratio, 
    medz.count.prop.ratio.top0, 
    medz.count.prop.ratio.top5)

count.ratio.peerless.plot = ggplot(data = medz.count.prop.ratio, aes(x = MAF, y = ESTIM)) +
    theme_bw() + xlab('Minor allele frequency (%)') + ylab('Enrichment') +
    geom_abline(intercept = 1, slope = 0) +
    geom_pointrange(aes(x = MAF, ymin = CI.LOW, ymax = CI.HIGH, colour = TYPE, alpha = PEER), position = position_dodge(width = 0.7)) +
    scale_colour_manual(values=colors.type) + guides(colour=FALSE) +
    scale_alpha_manual(values = peer.alphas, name = '') +
    facet_wrap( ~ TYPE, scales = "free") + 
    theme(strip.background = element_blank(), 
        legend.position = c(.2,.8), 
        legend.key = element_blank(), 
        legend.background = element_blank())
count.ratio.peerless.plot

# Save workspace image
save.image(paste0(dir, '/data/suppfig.count.enrichments.peer.effect.RData'))
