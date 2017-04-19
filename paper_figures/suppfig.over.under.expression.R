#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

## look into the expression distribution for overexpression and underexpression outliers

library(ggplot2)
library(cowplot)

## load figure panel with enrichments split by over and underexpression outliers
load(paste0(dir, '/data/figure3b.feature.enrichments.RData'))
load(paste0(dir, '/data/figure2c.ASE.enrichments.RData'))

fontsize = 11
fontsizes = theme(axis.text = element_text(size = fontsize),
                  axis.title = element_text(size = fontsize),
                  legend.text = element_text(size = fontsize - 1),
                  legend.title = element_text(size = fontsize - 1),
                  legend.key = element_blank(),
                  legend.background = element_blank(),
                  strip.text = element_text(face = 'plain', size = fontsize + 1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())

## read in median expression values
medians = read.table(paste0(dir, '/preprocessing/gtex_2015-01-12_median_rpkm.txt'), header = T, stringsAsFactors = F)
## read in outlier list
outliers.nothresh = read.table(paste0(dir, '/data/outliers_medz_nothreshold_picked.txt'), header = T, stringsAsFactors = F)

joined = merge(outliers.nothresh, medians, by.x = 'GENE', by.y = 'Gene')
joined$Type = ifelse(joined$Z > 0, 'Overexpression', 'Underexpression')
joined$Type[abs(joined$Z) < 2] = 'All'
extra = joined[joined$Type != 'All',]
extra$Type = 'All'
joined = rbind(joined, extra)

dist.plot = ggplot(joined, aes(log2(Median + 2), fill = Type)) + geom_density(alpha = 0.5, colour = 'white') +
    xlab(expression('Log'[2]~'(median RPKM + 2)')) + ylab('Density') +
    theme_bw() + scale_fill_manual(values = c('darkgrey','dodgerblue4','dodgerblue'), name = '') +
    fontsizes + theme(legend.position = c(0.7, 0.85))

add.effects.plot.over.under = add.effects.plot.over.under +
    fontsizes + theme(legend.position = c(0.2, 0.85))

ase.plot.rare = ase.plot.rare + fontsizes

## combine the two plots into a two-paneled figure
combined = ggdraw() +
    draw_plot(ase.plot.rare, 0.2, 0.5, 0.6, 0.5) +
    draw_plot(dist.plot, 0, 0, 5/11, 0.5) +
    draw_plot(add.effects.plot.over.under, 5/11, 0, 6/11, 0.5) +
    draw_plot_label(c('a', 'b', 'c'),
                    c(0.2, 0, 5/11),
                    c(1, 0.5, 0.5), size = 14)

pdf(paste0(dir, '/paper_figures/suppfig.over.under.expression.pdf'), height = 10, width = 11)

combined

dev.off()

## Wilcoxon tests comparing the log2(median RPKM+2) of genes with over(under)expression outliers to all genes
over = log2(joined$Median[joined$Type == 'Overexpression'] + 2)
under = log2(joined$Median[joined$Type == 'Underexpression'] + 2)
all = log2(joined$Median[joined$Type == 'All'] + 2)

print(wilcox.test(over, all))

print(wilcox.test(under, all))
