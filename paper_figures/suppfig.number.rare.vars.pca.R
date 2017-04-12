#!/usr/bin/env Rscript

## PATHS TO BE SET ###########
dir = Sys.getenv('RAREVARDIR')
#subjectFile = 'GTEx_Analysis_2015-01-12_Annotations_SubjectPhenotypesDS.txt' # adjust path as necessary
#covariatesDir = 'eQTLInputFiles/covariates/' # adjust path as necessary
subjectFile = '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/sample_annotations/GTEx_Analysis_2015-01-12_Annotations_SubjectPhenotypesDS.txt'
covariatesDir = '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_data/eQTLInputFiles/covariates/'


# looking into whether or not individuals with a lot of rare variants are outliers in PCA space

library(ggplot2)
library(reshape2)
library(cowplot)

# Generated when running figure2a.count.enrichments.R
load(paste0(dir, '/data/figure2a.count.enrichments.RData'))

# process variant counts (read in and stored in the RData file)
rownames(varCounts) = sub('-', '.', rownames(varCounts), fixed = T)
colnames(varCounts)[2] = "Indel"
varCounts.melted = melt(varCounts)

highRareCount = setdiff(rownames(varCounts), sub('-', '.', indIncl, fixed = T))

font.sizing = theme(axis.text = element_text(size = 12),
                    axis.title = element_text(size = 12),
                    legend.text = element_text(size = 12),
                    strip.background = element_blank(),
                    strip.text = element_text(size = 13),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank())

# histogram of the number of variants
options(scipen=999)
colours = c('dodgerblue3','springgreen4','orangered3')
counthist = ggplot(varCounts.melted, aes(x = value, fill = variable, colour = 'white')) + geom_histogram() +
    facet_grid(.~variable, scale = "free") + ylab('Number of individuals') + xlab('Number of variants') +
    theme_bw() + guides(fill = FALSE, colour = FALSE) +
    scale_fill_manual(values = colours) + scale_colour_manual(values = 'white') + font.sizing

# add vertical lines to indicate medians
medians = apply(varCounts, 2, median, na.rm = T)
medians.df = melt(medians)
medians.df$variable = factor(rownames(medians.df), levels = rownames(medians.df))
counthist = counthist + geom_vline(data = medians.df, aes(xintercept = value))

# update font sizing for enrichment plot
count.ratio.plot.subset = count.ratio.plot.subset + font.sizing

# read in covariates to extract genotype PCs
# need several tissue to cover all people (449 with genotype PCs available)
ms = read.table(paste0(covariatesDir, 'Muscle_Skeletal_Analysis.covariates.txt'), header = T, row.names = 1)
ms = t(ms[c(1:3), ])
wb = read.table(paste0(covariatesDir, 'Whole_Blood_Analysis.covariates.txt'), header = T, row.names = 1)
wb = t(wb[c(1:3), ])
cortex = read.table(paste0(covariatesDir, 'Brain_Frontal_Cortex_BA9_Analysis.covariates.txt'), header = T, row.names = 1)
cortex = t(cortex[c(1:3), ])
ad = read.table(paste0(covariatesDir, 'Adipose_Subcutaneous_Analysis.covariates.txt'), header = T, row.names = 1)
ad = t(ad[c(1:3), ])
th = read.table(paste0(covariatesDir, 'Thyroid_Analysis.covariates.txt'), header = T, row.names = 1)
th = t(th[c(1:3), ])
pcs = as.data.frame(unique(rbind(ms, wb, cortex, ad, th)))

# read in subject information (to extract ancestry info)
subj = read.delim(subjectFile, header = T, row.names = 1)[, c(4,5)]
# replace - by . to match the row names on pcs data frame
rownames(subj) = sub('-', '.', rownames(subj), fixed = T)
# replace codes by their meanings
subj$RACE[subj$RACE == 99] = 98
populations = c('White','White with WGS','Asian','Black','Native American','Unknown')
subj$RACE = factor(subj$RACE, levels = c(3,0,1,2,4,98), labels = populations)
subj$ETHNCTY[subj$ETHNCTY > 97] = 97
subj$ETHNCTY = factor(subj$ETHNCTY, levels = c(0,1,97), labels = c('Not Hispanic', 'Hispanic', 'Unknown'))

# add this race, ethnicity, and count info to PCs data frame
pcs$RACE = subj$RACE[match(rownames(pcs),rownames(subj))]
pcs$ETH = subj$ETHNCTY[match(rownames(pcs),rownames(subj))]
pcs$HIGHRARE = (rownames(pcs) %in% highRareCount)
# color by race, but make separate color for europeans with and without WGS
pcs$COL = pcs$RACE
pcs$COL[rownames(pcs) %in% rownames(varCounts)] = "White with WGS"
# order points for nicer plotting (with white with WGS plotted last)
pcs = pcs[order(as.character(pcs$COL)), ]

# plot PCs, color by race and circle for high number of variants
pca.colours = c('dodgerblue2','#64EBFF','darkolivegreen3','darkgoldenrod2','#FF3D3D','darkgrey')
names(pca.colours) = populations
coords = c(-0.025, -0.01, -0.05, 0.025)
boxdata = data.frame(x1 = c(rep(coords[1],3), rep(coords[2],2), coords[1]),
                 x2 = c(rep(coords[2],2), coords[1], coords[2], 0.14, 0.05),
                 y1 = c(coords[4], coords[3], rep(coords[4],3), coords[3]),
                 y2 = c(coords[4], rep(coords[3],3), -0.14, -0.59))
pca = ggplot(pcs, aes(x = C1, y = C2)) + geom_point(aes(colour = COL)) + theme_bw() +
    scale_colour_manual(values = pca.colours, name = element_blank()) + xlab('PC1') + ylab('PC2') +
    geom_point(data = pcs[pcs$HIGHRARE, ], aes(x = C1, y = C2), pch = 1, size = 2, stroke = 1, colour = "black") +
    geom_segment(data = boxdata, aes(x = x1, y = y1, xend = x2, yend = y2), linetype = 2, colour = 'darkgrey') +
    theme(legend.position = c(0.1, 0.5),
          legend.key = element_blank(),
          legend.background = element_blank()) + font.sizing

# plot zoom-in of top right corner of plot
pca.subset = ggplot(pcs, aes(x = C1, y = C2, colour = COL)) + geom_point() + theme_bw() +
    scale_colour_manual(values = pca.colours, name = element_blank()) + xlim(coords[1], coords[2]) + ylim(coords[3], coords[4]) + 
    geom_point(data = pcs[pcs$HIGHRARE, ], aes(x = C1, y = C2), pch = 1, size = 2, stroke = 1, color = "black") +
    guides(colour = FALSE) + xlab('') + ylab('') + font.sizing +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.margin = unit(c(1,1,1,-6), 'mm'))

# Plot of MAFs of of rare variants in individual european populations
# data generated by running ../feature_construction/assess.euro.subpop.R
load(paste0(dir, '/data/euro.subpop.RData'))

# combine SNVs and indels
euro.plotdata = rbind(indel.plotdata, snv.plotdata)

euro.maf = ggplot(euro.plotdata, aes(x = 100 * MAF, fill = Subpopulation)) +
    geom_histogram(binwidth = 0.1) +
    facet_grid(.~Subpopulation) + theme_bw() + guides(fill = FALSE) +
    ylab('Number of variants') + xlab("Minor allele frequency (%)") +
    scale_y_continuous(labels = scales::comma) + font.sizing


combined.plot = ggdraw() +
    draw_plot(counthist, 0.02, 12.1/14.5, 0.98, 2.4/14.5) +
    draw_plot(pca, 0, 5.5/14.5, 1, 6.6/14.5) +
    draw_plot(pca.subset, 0.48, 6.1/14.5, 0.47, 4.2/14.5) +
    draw_plot(count.ratio.plot.subset, 0, 2.5/14.5, 1, 3/14.5) +
    draw_plot(euro.maf, 0, 0, 1, 2.5/14.5) +
    draw_plot_label(c('a','b','c','d'), c(0,0,0,0), c(1, 12.1/14.5, 5.5/14.5, 2.5/14.5), size = 16)


pdf(paste0(dir, '/paper_figures/suppfig.varcounts.pca.pdf'), height = 14.5, width = 10.5)

print(combined.plot)

dev.off()
