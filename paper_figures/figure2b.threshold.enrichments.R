#!/usr/bin/env Rscript

# Load required packages
require(ggplot2)
require(reshape2)

dir = Sys.getenv('RAREVARDIR')

#----------------- FUNCTIONS

# for each variant MEDZ threshold,
# get the ratio of the portions of outliers/non outliers with rare variants
# and the associated confidence intervals (exp of Wald interval on the log of the proportion ratio)
proportion.ratios = function(counts, thresholds) {
    # make empty data frame and fill it from rbinds (not very many, so it's fine)
    results = data.frame(ESTIM=numeric(), CI.LOW=numeric(), CI.HIGH=numeric(), THRESH = numeric(), COUNT = numeric(),
        stringsAsFactors = F)
    for (thresh in thresholds) {
        props = proportion.ratios.helper(counts, thresh)
        if (!is.null(props)) {
            results = rbind(results, props)
        }
    }
    return(results)
}

# actually does the work of the main function for the given MEDZ or Z threshold
proportion.ratios.helper = function(counts, thresh, absval = TRUE) {
    if (absval) {
        counts.subset = counts[abs(counts$RANK) >= thresh, c('Y','any_variant')]
    } else {
        if (thresh <= 0) {
            counts.subset = counts[counts$RANK <= thresh, c('Y','any_variant')]
        } else {
            counts.subset = counts[counts$RANK >= thresh, c('Y','any_variant')]
        }
    }
    summary.counts = as.data.frame(table(counts.subset))
    # check that no zero values, and right number of values
    if (nrow(summary.counts) != 4 | min(summary.counts$Freq) == 0) {
        cat("Warning: Skipping this threshold because zeros in contingency table.\n")
        return(NULL)
    }
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
    dfrow = list(ESTIM=estimate, CI.LOW=min.ci, CI.HIGH=max.ci, THRESH=thresh, COUNT = sum(counts.subset$Y))
    return (dfrow)
}

read.counts = function(filename) {
    df = read.table(filename, header = T, stringsAsFactors = F)
    df$any_variant = (df$n_variants > 0) + 0
    return(df)
}

#----------------- MAIN
# Read in the count feature data for medz as well as for individual tissues (both with and without protein-coding variants)
medz = read.counts(paste0(dir, '/data/medz/10kb_features_SNPs_counts_MAF0-1_nothreshold.txt'))
tissuefiles = list.files(paste0(dir, '/data/singlez'), '10kb_features_SNPs_MAF0-1_nothreshold_*', full.names = T)
tissuenames = sub(paste0(dir, "/data/singlez/10kb_features_SNPs_MAF0-1_nothreshold_"), "", tissuefiles)
tissuenames = sub(".txt", "", tissuenames)
singlez = lapply(tissuefiles, read.counts)
names(singlez) = tissuenames

# Find the estimated proportions and CIs for varying MEDZ thresholds
thresholds = 1:5
medz.props = proportion.ratios(medz, thresholds)

# do the same for varying Z thresholds in single tissues
singlez.props = lapply(singlez, proportion.ratios, thresholds)

# combine the single tissues into a single data.frame
singlez.props.all = data.frame(ESTIM = numeric(), CI.LOW = numeric(), CI.HIGH = numeric(), THRESH = numeric(), COUNT = numeric(), TISSUE = character())
for (tis in names(singlez.props)) {
    tis.df = singlez.props[[tis]]
    tis.df$TISSUE = tis
    singlez.props.all = rbind(singlez.props.all, tis.df)
}

# also try combining the single tissues into one data frame *before* calculating the enrichments
thresholds.single = 1:10
singlez.all = data.frame(gene_id = character(), n_variants = numeric(), ID = character(), RANK = numeric(), Y = numeric(), any_variant = numeric())
for (tis in names(singlez)) {
    singlez.all = rbind(singlez.all, singlez[[tis]])
}
singlez.combined.props = proportion.ratios(singlez.all, thresholds.single)

# merge the combined single tissue with the multi-tissue results
singlez.combined.props$type = "Single-tissue outlier"
medz.props$type = "Multi-tissue outlier"
both.props = rbind(singlez.combined.props, medz.props)

## plot for main text with combined single tissue and also multi-tissue outliers
colors = c('dodgerblue3', 'mediumorchid4')
names(colors) = c("Multi-tissue outlier","Single-tissue outlier")

both.props$labely = ifelse(both.props$type == "Multi-tissue outlier", yes = both.props$CI.HIGH + 0.05, no = 0.95)
index.above = both.props$type == "Single-tissue outlier" & both.props$THRESH %in% c(2,4,6,8)
both.props$labely[index.above] = both.props$CI.HIGH[index.above] + 0.07

props.plot.medz.singlez = ggplot(both.props, aes(x = THRESH, y = ESTIM, colour = type)) +
    geom_pointrange(aes(ymin = CI.LOW, ymax = CI.HIGH)) + theme_bw() +
    xlab('Z-score threshold') + ylab('Enrichment') + geom_hline(yintercept = 1) +
    scale_colour_manual(values = colors, name = "") +
    theme(legend.key = element_blank(),
          legend.position = c(0.8,0.65),
          legend.background = element_blank()) +
    ylim(c(min(both.props$CI.LOW) - 0.1, max(both.props$CI.HIGH) + 0.1)) +
    annotate('text', x = both.props$THRESH, y = both.props$labely,
             label = formatC(both.props$COUNT, format = "d", big.mark = ","), size = 3.5) +
    scale_x_continuous(breaks = thresholds.single,
                       labels = thresholds.single,
                       limits = c(thresholds.single[1] - .25, thresholds.single[length(thresholds.single)] + 0.25))

## Supplementary figure
## plots for individual single tissues
# read in file with tissue names in long form
tissues = read.table('gtex_tissue_colors.txt', header = T, stringsAsFactors = F, sep = "\t")
tissues$names_pretty = sub(" - ", " -\n", tissues$tissue_site_detail, fixed = T)
tissues$names_pretty = sub(" \\(.*?\\)", "", tissues$names_pretty)
tissues$names_pretty[tissues$names_pretty == "Heart -\nAtrial Appendage"] = "Heart - Atrial\nAppendage"
tissues$names_pretty[tissues$names_pretty == "Skin -\nNot Sun Exposed"] = "Skin - Not\nSun Exposed"
tissues$names_pretty[tissues$names_pretty == "Brain -\nNucleus accumbens"] = "Brain - Nucleus\naccumbens"
tissues$names_pretty[tissues$names_pretty == "Brain -\nCerebellar Hemisphere"] = "Brain - Cerebellar\nHemisphere"
tissues$names_pretty[tissues$names_pretty == "Brain -\nAnterior cingulate cortex"] = "Brain - Anterior\ncingulate cortex"
tissues$names_pretty[tissues$names_pretty == "Cells -\nEBV-transformed lymphocytes"] = "Cells -\nLymphocytes"
tissues$names_pretty[tissues$names_pretty == "Cells -\nTransformed fibroblasts"] = "Cells -\nFibroblasts"
tissues$names_pretty[tissues$names_pretty == "Esophagus -\nGastroesophageal Junction"] = "Esophagus -\nJunction"
singlez.props.all$TISSUE_NAME = tissues$names_pretty[match(singlez.props.all$TISSUE, tissues$tissue_site_detail_id)]

single.tissues.plot = ggplot(singlez.props.all, aes(x = THRESH, y = ESTIM)) + geom_hline(yintercept = 1) +
    geom_pointrange(aes(ymin = CI.LOW, ymax = CI.HIGH), colour = "mediumorchid4", size = 0.3) +
    theme_bw() + xlab(expression("|Z-score|")) + ylab("Enrichment") +
    facet_wrap(~TISSUE_NAME, nrow = 7, scales = "free_y") +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 8),
          panel.spacing.x = unit(0.8, 'pt'),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

# Save workspace image
save.image(paste0(dir, '/data/figure2b.threshold.enrichments.RData'))


