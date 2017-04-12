#!/usr/bin/env Rscript

dir = Sys.getenv('RAREVARDIR')

## Supplemental figure to accompany figure 2
## Includes threshold enrichemnts for individual tissues and
## count enrichments by MAF with and without protein-coding regions.

library(ggplot2)
library(cowplot)

## load figure panel with enrichments split by over and underexpression outliers
load(paste0(dir, '/data/figure2a.count.enrichments.RData'))
load(paste0(dir, '/data/figure2b.threshold.enrichments.RData'))

fontsize = 8
fontsizes = theme(axis.text = element_text(size = fontsize),
                  axis.title = element_text(size = fontsize + 1),
                  legend.text = element_text(size = fontsize),
                  legend.title = element_text(size = fontsize),
                  legend.key = element_blank(),
                  legend.background = element_blank(),
                  strip.text = element_text(size = fontsize + 1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())

count.ratio.plot.with.withoutPC = count.ratio.plot.with.withoutPC + fontsizes
single.tissues.plot = single.tissues.plot + fontsizes

## combine the two plots into a two-paneled figure
combined = ggdraw() +
    draw_plot(count.ratio.plot.with.withoutPC, 0.1, 4/5, 0.8, 1/5) +
    draw_plot(single.tissues.plot, 0, 0, 1, 4/5) +
    draw_plot_label(c('a', 'b'), c(0.1, 0), c(1, 4/5), size = 11)

pdf(paste0(dir, '/paper_figures/suppfig.count.enrichments.pdf'), height = 11.5, width = 10)

combined

dev.off()

