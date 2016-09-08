#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

# Load required packages and functions
require(ggplot2)
require(gtable)
source('enrichment.functions.R')
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

# Read in features and counts
MAFs = c('0-1', '1-5', '5-10', '10-25')
TYPEs = c('SNPs', 'indels', 'HallLabSV')
type.names = c('SNV','Indel','SV')

OUT.DIR = paste0(dir, '/data/medz/')
medz.count.withPC = read_file(OUT.DIR, MAFs, TYPEs, '_counts_MAF', type.names, prefix = '10kb_')
medz.count.withoutPC = read_file(OUT.DIR, MAFs, TYPEs, '_counts_MAF', type.names, prefix = '10kb_noPC_')

# Turn counts into binary presence/absence 
medz.count.withPC$any_variant = (medz.count.withPC$n_variants > 0) + 0
medz.count.withoutPC$any_variant = (medz.count.withoutPC$n_variants > 0) + 0

# Make plot of count enrichments by MAF faceted by variant type
# Include error bars
colors.type = c('dodgerblue3','springgreen4','orangered3')
names(colors.type) = type.names

# Make plot of count enrichment for SNPs only using ratio of proporions instead of log odds
medz.count.prop.ratio.withPC = proportion.ratios(medz.count.withPC)
medz.count.prop.ratio.withPC$MAF = factor(medz.count.prop.ratio.withPC$MAF, levels = MAFs)
medz.count.prop.ratio.withPC$TYPE = factor(medz.count.prop.ratio.withPC$TYPE, levels = type.names)

medz.count.prop.ratio.withoutPC = proportion.ratios(medz.count.withoutPC)
medz.count.prop.ratio.withoutPC$MAF = factor(medz.count.prop.ratio.withoutPC$MAF, levels = MAFs)
medz.count.prop.ratio.withoutPC$TYPE = factor(medz.count.prop.ratio.withoutPC$TYPE, levels = type.names)

count.ratio.plot.withPC = ggplot(data = medz.count.prop.ratio.withPC, aes(x = MAF, y = ESTIM)) +
    theme_bw() + xlab('Minor allele frequency (%)') + ylab('Enrichment') +
    geom_abline(intercept = 1, slope = 0) +
    geom_pointrange(aes(x = MAF, ymin = CI.LOW, ymax = CI.HIGH, colour = TYPE)) +
    scale_colour_manual(values=colors.type) + guides(colour=FALSE) +
    facet_wrap( ~ TYPE, scales = "free") + 
    theme(strip.background = element_blank())

# make same plot (with PC) but excluding individuals with more than 60,000 variants
# for supplemental figure
varCounts = read.table(paste0(dir, '/data/variant_counts_per_individual_all.txt'), header = T, row.names = 1)
indIncl = rownames(varCounts)[rowSums(varCounts, na.rm = T) <= 60000]
medz.count.withPC.subset = medz.count.withPC[medz.count.withPC$ID %in% indIncl,]
medz.count.prop.ratio.withPC.subset = proportion.ratios(medz.count.withPC.subset)
medz.count.prop.ratio.withPC.subset$MAF = factor(medz.count.prop.ratio.withPC.subset$MAF, levels = MAFs)
medz.count.prop.ratio.withPC.subset$TYPE = factor(medz.count.prop.ratio.withPC.subset$TYPE, levels = type.names)
medz.count.prop.ratio.withPC.subset$GROUP = "Excluding individuals\nwith > 60000 variants"
medz.count.prop.ratio.withPC$GROUP = "All variants"
medz.count.prop.ratio.withPC.subset = rbind(medz.count.prop.ratio.withPC.subset, medz.count.prop.ratio.withPC)

alphas = c(1, 0.4)
count.ratio.plot.withPC.subset = ggplot(data = medz.count.prop.ratio.withPC.subset,
                                        aes(x = MAF, y = ESTIM, alpha = GROUP, colour = TYPE)) +
    theme_bw() + xlab('Minor allele frequency (%)') + ylab('Enrichment') +
    geom_abline(intercept = 1, slope = 0) +
    geom_pointrange(aes(x = MAF, ymin = CI.LOW, ymax = CI.HIGH, group = interaction(MAF, GROUP)),
                    position = position_dodge(width = 0.4)) +
    scale_colour_manual(values = colors.type) + guides(colour = FALSE) +
    scale_alpha_manual(values = alphas, name = NULL) +
    facet_wrap( ~ TYPE, scales = "free") +
    theme(strip.background = element_blank(),
          legend.key = element_blank())

# combine with and without PC in the same plot
medz.count.prop.ratio.withoutPC$GROUP = "Excluding\nexon variants"
medz.count.prop.ratio.all = rbind(medz.count.prop.ratio.withPC, medz.count.prop.ratio.withoutPC)
medz.count.prop.ratio.all$GROUP = factor(medz.count.prop.ratio.all$GROUP)

alphas = c(1, 0.4)
fsize = 9
count.ratio.plot.with.withoutPC = ggplot(medz.count.prop.ratio.all, aes(x = MAF, y = ESTIM, alpha = GROUP)) +
    theme_bw() + xlab('Minor allele frequency (%)') + ylab('Enrichment') +
    geom_abline(intercept = 1, slope = 0) +
    geom_pointrange(aes(x = MAF, ymin = CI.LOW, ymax = CI.HIGH, colour = TYPE)) +
    scale_colour_manual(values = colors.type) + guides(colour = FALSE) +
    scale_alpha_manual(values = alphas, name = NULL) +
    facet_wrap( ~ TYPE, scales = "free") +
    theme(strip.text.x = element_text(size = fsize + 1),
          strip.background = element_blank(),
          legend.key = element_blank(),
          axis.text = element_text(size = fsize),
          axis.title = element_text(size = fsize),
          legend.text = element_text(size = fsize),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

# output the supplementary figure
pdf(paste0(dir, '/paper_figures/suppfig.figure2a.withoutPC.pdf'), height = 2.3, width = 7.2)

count.ratio.plot.with.withoutPC

dev.off()

# Save workspace image
save.image(paste0(dir, '/data/figure2a.count.enrichments.RData'))

